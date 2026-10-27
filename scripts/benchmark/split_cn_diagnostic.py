#!/usr/bin/env python
"""Split-leakage diagnostic (benchmark-label-free; DESIGN round-3 #16): how similar are TEST regulators' CN
profiles to TRAIN regulators' across the DepMap S2 panel?

Uses only the split (results/split/), the universe regulator lists and DepMap CN. No benchmark label is read.
x = log2(clip(relCN, 2^-2, 2^3)) (DESIGN §3); X genes and Xp/Xq arms are sex-centred (median within
Model sex; Unknown/missing pooled), as for the frozen variable-CN flag.

Per line and TEST regulator: max |Pearson r| with any TRAIN regulator, the arm of that TRAIN regulator, and
whether it is the TEST arm's sibling (other arm of the same chromosome). The reference is the same statistic
for TRAIN fold-2 regulators against TRAIN fold-1 regulators (the 2-fold tuning split also splits siblings).
Arm level: |r| between each TEST arm's arm CN and every TRAIN arm.

Outputs: results/split/cn_leakage_regulators_<line>.tsv, results/split/cn_leakage_arms.tsv,
results/split/cn_leakage_summary.json
Usage: .venv/bin/python scripts/benchmark/split_cn_diagnostic.py
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from depmaptr import benchmark as bm

ROOT = bm.ROOT
PROC = ROOT / "data" / "processed"
OUT = ROOT / "results" / "split"
THRESH = (0.5, 0.8, 0.9, 0.95)


def sibling(arm: str) -> str:
    return arm[:-1] + ("q" if arm[-1] == "p" else "p")


def standardise(x: np.ndarray) -> np.ndarray:
    x = x - x.mean(0)
    sd = x.std(0)
    sd[sd == 0] = np.inf
    return x / sd


def sex_centre(x: pd.DataFrame, sex: np.ndarray) -> pd.DataFrame:
    x = x.copy()
    for s in np.unique(sex):
        rows = sex == s
        x.loc[rows] = x.loc[rows] - x.loc[rows].median()
    return x


def max_corr(za: np.ndarray, zb: np.ndarray, chunk: int = 1024):
    n = za.shape[0]
    best = np.empty(za.shape[1])
    arg = np.empty(za.shape[1], dtype=np.int64)
    for i in range(0, za.shape[1], chunk):
        with np.errstate(all="ignore"):                  # spurious Accelerate flags (DESIGN §10)
            c = np.abs(za[:, i:i + chunk].T @ zb) / n
        arg[i:i + chunk] = c.argmax(1)
        best[i:i + chunk] = c.max(1)
    return best, arg


def summarise(v: np.ndarray) -> dict:
    return {"n": int(len(v)), "median": float(np.median(v)), "q90": float(np.quantile(v, 0.9)),
            **{f"frac_gt_{t}": float(np.mean(v > t)) for t in THRESH}}


def main() -> int:
    panel = pd.read_csv(PROC / "panel_s2.tsv", sep="\t")
    sex = panel["sex"].fillna("Unknown").to_numpy()
    units = pd.read_csv(OUT / "arm_units.tsv", sep="\t")
    asplit = dict(zip(units["arm"], units["split"]))
    summary = {"definition": __doc__.split("\n\n")[1], "lines": {}}
    for line in bm.LINES:
        regs = pd.read_csv(bm.bench_dir(ROOT) / f"universe_regulators_{line}.tsv", sep="\t", dtype={"chrom": str})
        cn = pd.read_parquet(PROC / "depmap_cn_matched.parquet", columns=["ModelID", *regs["regulator"]])
        cn = cn.set_index("ModelID").loc[panel["ModelID"]]
        x = np.log2(cn.clip(2.0 ** -2, 2.0 ** 3))
        x = x.fillna(x.median())
        isx = (regs["chrom"] == "X").values
        x.loc[:, regs["regulator"][isx].values] = sex_centre(x.loc[:, regs["regulator"][isx].values], sex)
        z = standardise(x.to_numpy(np.float64))
        te = np.flatnonzero((regs["split"] == "TEST").values)
        tr = np.flatnonzero((regs["split"] == "TRAIN").values)
        f1 = np.flatnonzero(((regs["split"] == "TRAIN") & (regs["fold"] == 1)).values)
        f2 = np.flatnonzero(((regs["split"] == "TRAIN") & (regs["fold"] == 2)).values)
        best, arg = max_corr(z[:, te], z[:, tr])
        rows = regs.iloc[te][["regulator", "arm"]].copy()
        rows["max_abs_r_train"] = best
        rows["nearest_train_arm"] = regs["arm"].values[tr][arg]
        rows["nearest_is_sibling"] = rows["nearest_train_arm"] == rows["arm"].map(sibling)
        rows["sibling_in_train"] = rows["arm"].map(sibling).map(asplit).eq("TRAIN")
        rows.to_csv(OUT / f"cn_leakage_regulators_{line}.tsv", sep="\t", index=False, float_format="%.4f")
        bref, _ = max_corr(z[:, f2], z[:, f1])
        summary["lines"][line] = {
            "test_vs_train": summarise(best),
            "test_vs_train_sibling_in_train": summarise(best[rows["sibling_in_train"].values]),
            "test_vs_train_sibling_not_in_train": summarise(best[~rows["sibling_in_train"].values]),
            "frac_nearest_is_sibling": float(rows["nearest_is_sibling"].mean()),
            "reference_fold2_vs_fold1": summarise(bref)}
        print(line, json.dumps(summary["lines"][line]["test_vs_train"]), file=sys.stderr)

    arm = pd.read_parquet(PROC / "depmap_arm_cn_derived.parquet")
    a = arm.pivot(index="ModelID", columns="arm_id", values="median_log2_cn").loc[panel["ModelID"]]
    a = a[[c for c in units["arm"] if c in a.columns]]
    xa = [c for c in a.columns if c.startswith("X")]
    a.loc[:, xa] = sex_centre(a.loc[:, xa], sex)
    a = a.fillna(a.median())
    za = standardise(a.to_numpy(np.float64))
    cols = list(a.columns)
    with np.errstate(all="ignore"):
        C = np.abs(za.T @ za) / za.shape[0]
    rows = []
    for i, arm_i in enumerate(cols):
        if asplit.get(arm_i) != "TEST":
            continue
        trj = [j for j, c in enumerate(cols) if asplit.get(c) == "TRAIN"]
        j = trj[int(np.argmax(C[i, trj]))]
        sib = sibling(arm_i)
        rows.append({"test_arm": arm_i, "sibling": sib, "sibling_split": asplit.get(sib, "absent"),
                     "abs_r_sibling": float(C[i, cols.index(sib)]) if sib in cols else np.nan,
                     "max_abs_r_train_arm": float(C[i, j]), "nearest_train_arm": cols[j]})
    arms = pd.DataFrame(rows)
    arms.to_csv(OUT / "cn_leakage_arms.tsv", sep="\t", index=False, float_format="%.4f")
    summary["arms"] = {"n_test_arms_with_train_sibling": int((arms["sibling_split"] == "TRAIN").sum()),
                       "chromosomes_split_across": sorted({r[:-1] for r in arms.loc[arms["sibling_split"] == "TRAIN",
                                                                                     "test_arm"]}, key=bm.chrom_order),
                       "median_abs_r_with_train_sibling": float(arms.loc[arms["sibling_split"] == "TRAIN",
                                                                         "abs_r_sibling"].median())}
    (OUT / "cn_leakage_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary["arms"], indent=1))
    print(arms.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
