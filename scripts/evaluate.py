#!/usr/bin/env python
"""Evaluate score files on the TRAIN benchmark (DESIGN §6, §8.2-8.3). Every call is appended to
results/tuning_log.tsv by depmaptr.benchmark.evaluate(). TEST is refused here; only scripts/final_eval.py
can evaluate it.

Examples
  .venv/bin/python scripts/evaluate.py M1                       # HCT116 TRAIN, all metrics
  .venv/bin/python scripts/evaluate.py M1 B3 --line both
  .venv/bin/python scripts/evaluate.py M1_g1 --compare M1      # paired comparison + §8.2 adoption rule
  .venv/bin/python scripts/evaluate.py M1 --subset class_i      # secondary universe: own-gene class (i)
RUN is a run id under data/results/scores/ or a directory holding scores.npz.
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from depmaptr import benchmark as bm

SUBSETS = {"class_i": ("own_class", "i"), "class_ii": ("own_class", "ii"), "class_iii": ("own_class", "iii")}


def headline(out: dict) -> dict:
    m = out["metrics"]
    h = {"run_id": out["run_id"], "line": out["line"], "split": out["split"], "subset": out.get("subset"),
         "n_pairs": m["n_pairs"], "K": m[f"top{bm.PRIMARY_K:g}"]["K"], "OR@0.1%": m["primary_or"],
         "OR@0.01%": m["top0.0001"]["or"], "OR@1%": m["top0.01"]["or"],
         "fisher_p(indep.)": m[f"top{bm.PRIMARY_K:g}"]["fisher_p_greater"], "folds_log_or": m.get("folds_log_or"),
         "jackknife_ci95_OR": [float(np.exp(x)) for x in m["jackknife_log_or"].get("ci95", [np.nan, np.nan])]}
    if "perm_regulator" in m:
        h.update({"perm_p": m["perm_regulator"]["perm_p"], "perm_p_degree": m["perm_degree"]["perm_p"],
                  "OR/median_null": m["perm_regulator"]["or_over_median_null"],
                  "CMH_OR(reg x expr)": m["cmh_regulator_x_exprdecile"]["cmh_or"],
                  "sign_concordant_OR": m["sign"]["sign_concordant_or"],
                  "concordance": m["sign"].get("concordance"), "c0": m["sign"].get("c0"),
                  "sign_perm_p": m["sign"].get("sign_perm_p")})
    h["eval_wall_s"] = round(out["wall_s"], 2)
    return h


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--line", default=bm.PRIMARY_LINE, choices=list(bm.LINES) + ["both"])
    ap.add_argument("--split", default="TRAIN", help="only TRAIN is accepted here")
    ap.add_argument("--n-perm", type=int, default=bm.N_PERM)
    ap.add_argument("--n-boot", type=int, default=bm.N_BOOT)
    ap.add_argument("--quick", action="store_true", help="ORs, folds and jackknife only (no permutations)")
    ap.add_argument("--compare", metavar="INCUMBENT", help="paired comparison of each RUN against INCUMBENT")
    ap.add_argument("--subset", choices=sorted(SUBSETS))
    ap.add_argument("--variant", help="variant id for the log (default: from meta.json)")
    a = ap.parse_args()
    if a.split != "TRAIN":
        raise SystemExit("evaluate.py evaluates TRAIN only; TEST is evaluated once by scripts/final_eval.py")
    lines = list(bm.LINES) if a.line == "both" else [a.line]
    for line in lines:
        for run in a.runs:
            if a.compare:
                print(json.dumps(bm.compare(run, a.compare, line, "TRAIN", n_boot=a.n_boot), indent=1))
                continue
            mask = None
            if a.subset:
                col, val = SUBSETS[a.subset]
                regs = bm.load_bench(line, "TRAIN").regs
                mask = (regs[col].astype(str) == val).values
            out = bm.evaluate(run, line, "TRAIN", variant=a.variant, n_perm=a.n_perm, n_boot=a.n_boot,
                              full=not a.quick, reg_subset=mask, subset_name=a.subset)
            print(json.dumps(headline(out), default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
