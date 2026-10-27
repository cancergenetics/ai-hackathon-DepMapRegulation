#!/usr/bin/env python
"""Build the frozen benchmark tables and the held-out split (DESIGN v3.1 §6; docs/INTERFACES.md
"Benchmark tables"; definitions and implementation choices in docs/BENCHMARK.md).

Inputs (all label-free except the PIE caches):
  data/processed/pie_gene_map.tsv, gene_annotation.tsv, frozen_flags.tsv, panel_s2.tsv,
  depmap_expr_matched.parquet, starting_data/Model.csv (context-model assertion),
  data/processed/pie/{targets,features}_<line>.tsv, vocab_*.tsv, tested_<line>.parquet (own-gene tested
  status), loose_<line>.parquet (every PIE row with p_adj < 0.1 at full precision: positives, own-gene
  knockdown).

Outputs:
  data/processed/bench/universe_regulators_<line>.tsv   (TEST rows: own_class / degree_decile left empty)
  data/processed/bench/universe_targets_<line>.tsv
  data/processed/bench/positives_<line>_{TRAIN,TEST}.parquet
  data/processed/bench/regulator_labels_<line>_TEST.tsv (label-derived TEST columns; sealed with TEST)
  data/processed/bench/manifest.json
  results/split/regulator_split.tsv, arm_units.tsv, train_folds.tsv, benchmark_summary.json

Never run after the freeze: the TEST label files are sealed by scripts/benchmark/seal_test.sh.
Usage: .venv/bin/python scripts/benchmark/build_benchmark.py [--force]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

from depmaptr import benchmark as bm

ROOT = bm.ROOT
PROC = ROOT / "data" / "processed"
PIE = PROC / "pie"
OUT = bm.bench_dir(ROOT)
SPLIT_DIR = ROOT / "results" / "split"
# DESIGN §5.3 names "HCT 116 / HEK 293"; Model.csv CellLineName spells HEK-293, so the stripped name is used
CONTEXT_NAMES = {"ACH-000971": "HCT116", "ACH-001085": "HEK293"}


def cytoband_fill(ann: pd.DataFrame) -> pd.DataFrame:
    """Gene start/end; genes placed only by an HGNC band (no Ensembl coordinates) get that band's UCSC hg38
    extent (used only for the 5 Mb target bin; docs/BENCHMARK.md)."""
    pos = ann.set_index("depmap_column")[["start", "end", "chrom", "hgnc_location"]].copy()
    cb = pd.read_csv(ROOT / "data" / "external" / "annotation" / "cytoBand.hg38.txt.gz", sep="\t", header=None,
                     names=["chrom", "start", "end", "band", "stain"])
    cb["loc"] = cb["chrom"].str.removeprefix("chr") + cb["band"]
    ext = cb.set_index("loc")[["start", "end"]]
    miss = pos["start"].isna() & pos["hgnc_location"].notna()
    for g in pos.index[miss]:
        loc = str(pos.at[g, "hgnc_location"])
        if loc in ext.index:
            pos.at[g, "start"], pos.at[g, "end"] = ext.at[loc, "start"], ext.at[loc, "end"]
    return pos[["start", "end"]]


def own_gene_status(line: str, regs: pd.DataFrame, loose: pd.DataFrame) -> pd.DataFrame:
    """Own feature (same PIE id as the perturbed gene) tested? knocked down (p_adj < 0.05, LFC < 0)?"""
    vt = pd.read_csv(PIE / "vocab_targets.tsv", sep="\t")
    vf = pd.read_csv(PIE / "vocab_features.tsv", sep="\t")
    tcode = dict(zip(vt["target"], vt["code"]))
    fcode = dict(zip(vf["feature"], vf["code"]))
    pairs = pd.DataFrame({"pie_id": regs["pie_id"]})
    pairs["t"] = pairs["pie_id"].map(tcode)
    pairs["f"] = pairs["pie_id"].map(fcode)
    if pairs["t"].isna().any():
        raise SystemExit("universe regulator missing from the PIE target vocabulary")
    q = pairs.dropna(subset=["f"]).astype({"t": "int64", "f": "int64"})
    own = pl.from_pandas(q[["t", "f"]]).with_columns(pl.col("t").cast(pl.UInt16), pl.col("f").cast(pl.UInt16))
    hit = (pl.scan_parquet(PIE / f"tested_{line}.parquet").select(["t", "f", "lfc"])
           .join(own.lazy(), on=["t", "f"], how="inner").collect(engine="streaming").to_pandas())
    tested = pairs.merge(hit, on=["t", "f"], how="left")
    tested["own_tested"] = tested["lfc"].notna()
    lo = loose[loose["target"] == loose["feature"]][["target", "p_adj", "log2_fold_change"]]
    tested = tested.merge(lo, left_on="pie_id", right_on="target", how="left")
    tested.loc[tested["own_tested"] & tested["p_adj"].isna(), "p_adj"] = 1.0       # p_adj >= 0.1
    tested["own_class"] = bm.own_gene_class(tested["own_tested"], tested["p_adj"], tested["lfc"])
    return tested.set_index("pie_id")[["own_tested", "lfc", "p_adj", "own_class"]]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="overwrite existing benchmark tables")
    args = ap.parse_args()
    if (OUT / "manifest.json").exists() and not args.force:
        raise SystemExit(f"{OUT} already built; use --force (never after the freeze)")
    OUT.mkdir(parents=True, exist_ok=True)
    SPLIT_DIR.mkdir(parents=True, exist_ok=True)
    inputs = [PROC / "pie_gene_map.tsv", PROC / "gene_annotation.tsv", PROC / "frozen_flags.tsv",
              PROC / "panel_s2.tsv", PROC / "depmap_expr_matched.parquet"]
    for line in bm.LINES:
        inputs += [PIE / f"targets_{line}.tsv", PIE / f"features_{line}.tsv", PIE / f"loose_{line}.parquet",
                   PIE / f"tested_{line}.parquet"]
    inputs += [PIE / "vocab_targets.tsv", PIE / "vocab_features.tsv",
               ROOT / "data" / "external" / "annotation" / "cytoBand.hg38.txt.gz"]
    manifest = {"inputs": {str(p.relative_to(ROOT)): bm.sha256_file(p) for p in inputs},
                "seed": bm.SEED, "lines": {}}

    pie_map = pd.read_csv(PROC / "pie_gene_map.tsv", sep="\t", low_memory=False, dtype={"chrom": str})
    flags = pd.read_csv(PROC / "frozen_flags.tsv", sep="\t", dtype={"chrom": str})
    flags = bm.normalise_flags(flags)
    ann = pd.read_csv(PROC / "gene_annotation.tsv", sep="\t", low_memory=False, dtype={"chrom": str})
    gene_pos = cytoband_fill(ann)
    panel = pd.read_csv(PROC / "panel_s2.tsv", sep="\t")

    # ---- invariants (DESIGN §8.6, benchmark-relevant) ----
    sym = flags.reset_index().set_index("symbol")
    for g, arm in (("ERBB2", "17q"), ("MYC", "8q"), ("CDKN2A", "9p")):
        if sym.loc[g, "arm_id"] != arm:
            raise SystemExit(f"invariant: {g} not on {arm}")
    bm.pie_key_map(pie_map)                                  # asserts 1:1 PIE id <-> DepMap key
    models = pd.read_csv(ROOT / "starting_data" / "Model.csv", usecols=["ModelID", "StrippedCellLineName"])
    for mid, name in CONTEXT_NAMES.items():
        got = models.loc[models["ModelID"] == mid, "StrippedCellLineName"].tolist()
        if got != [name]:
            raise SystemExit(f"context model {mid}: Model.csv StrippedCellLineName {got} != {name}")
        if mid not in set(panel["ModelID"]):
            raise SystemExit(f"context model {mid} not in panel S2")

    # ---- DepMap expression (label-free target strata, regulator context expression) ----
    ex = pd.read_parquet(PROC / "depmap_expr_matched.parquet").set_index("ModelID")
    ex_s2 = ex.loc[panel["ModelID"]]
    expr_mean = ex_s2.mean(axis=0)
    expr_ctx = {line: ex.loc[bm.CONTEXT_MODEL[line]] for line in bm.LINES}

    regs_all, tgts_all, pos_all = {}, {}, {}
    for line in bm.LINES:
        lt = pd.read_csv(PIE / f"targets_{line}.tsv", sep="\t")
        ft = pd.read_csv(PIE / f"features_{line}.tsv", sep="\t")
        cons = ft[ft["ref_mean"] >= bm.CONSISTENT_REF_MEAN]
        if not (cons["n_targets"] == len(lt)).all():
            raise SystemExit(f"{line}: a consistently tested feature is not tested for every perturbation")
        regs, rc = bm.build_universe_regulators(pie_map, flags.reset_index(), lt["target"])
        tgts, tc = bm.build_universe_targets(pie_map, flags.reset_index(), cons["feature"], gene_pos)
        loose = pd.read_parquet(PIE / f"loose_{line}.parquet")
        own = own_gene_status(line, regs, loose)
        med_lfc = float(own.loc[own["own_tested"], "lfc"].median())
        if not med_lfc < 0:
            raise SystemExit(f"invariant: median on-target LFC in {line} is {med_lfc} (expected < 0)")
        regs["own_class"] = regs["pie_id"].map(own["own_class"])
        regs["own_lfc"] = regs["pie_id"].map(own["lfc"])           # build-time only; not written out
        regs["cn_stat"] = regs["regulator"].map(flags["n_lines_abs_log2cn_ge1"])
        if regs["cn_stat"].isna().any():
            raise SystemExit("universe regulator without a CN-variability statistic")
        regs["cn_tercile"] = bm.quantile_bins(regs["cn_stat"].values, 3)
        regs["expr_context"] = regs["regulator"].map(expr_ctx[line])        # NaN: CN-only regulator
        tgts["expr_mean_s2"] = tgts["target"].map(expr_mean)
        if tgts["expr_mean_s2"].isna().any():
            raise SystemExit("universe target without DepMap expression")
        tgts["expr_decile"] = bm.quantile_bins(tgts["expr_mean_s2"].values, 10)
        tgts["expr_hct116"] = tgts["target"].map(expr_ctx["hct116"])
        tgts["expr_hek293"] = tgts["target"].map(expr_ctx["hek293t"])
        pos = bm.build_positives(loose[loose["target"].isin(regs["pie_id"])], regs, tgts)
        # every positive is a trans universe pair
        rch = dict(zip(regs["regulator"], regs["chrom"]))
        tch = dict(zip(tgts["target"], tgts["chrom"]))
        assert (pos["regulator"].map(rch) != pos["target"].map(tch)).all()
        regs_all[line], tgts_all[line], pos_all[line] = regs, tgts, pos
        n_pairs = int(sum(len(tgts) - (tgts["chrom"] == c).sum() for c in regs["chrom"]))
        manifest["lines"][line] = {"regulator_counts": rc, "target_counts": tc, "n_pairs": n_pairs,
                                   "n_cn_only_regulators": int(regs["expr_context"].isna().sum())}
        print(f"{line}: {len(regs)} regulators, {len(tgts)} targets, {n_pairs} pairs, "
              f"{int(pos['pos_primary'].sum())} positives", file=sys.stderr)

    # ---- split (arm triplets; DESIGN §6) ----
    arms = pd.DataFrame({"arm": sorted(set(regs_all["hct116"]["arm"]) | set(regs_all["hek293t"]["arm"]),
                                       key=bm.arm_sort_key)})
    for line in bm.LINES:
        arms[f"n_{line}"] = arms["arm"].map(regs_all[line]["arm"].value_counts()).fillna(0).astype(int)
    units = bm.make_split(arms)
    units.to_csv(SPLIT_DIR / "arm_units.tsv", sep="\t", index=False)
    units.loc[units["split"] == "TRAIN", ["arm", "order", "fold", "n_hct116", "n_hek293t"]].to_csv(
        SPLIT_DIR / "train_folds.tsv", sep="\t", index=False)
    asplit = dict(zip(units["arm"], units["split"]))
    afold = dict(zip(units["arm"], units["fold"]))
    genes = pd.concat([regs_all[l][["regulator", "arm", "chrom"]] for l in bm.LINES]).drop_duplicates()
    genes = genes.sort_values("regulator").reset_index(drop=True)
    for line in bm.LINES:
        genes[f"in_{line}"] = genes["regulator"].isin(regs_all[line]["regulator"])
    genes["split"] = genes["arm"].map(asplit)
    genes["fold"] = genes["arm"].map(afold)
    if genes["split"].isna().any():
        raise SystemExit("regulator on an arm without a split")
    genes.rename(columns={"regulator": "gene"}).to_csv(SPLIT_DIR / "regulator_split.tsv", sep="\t", index=False)

    summary = {"n_arm_units": len(units), "arms_by_split": units["split"].value_counts().to_dict(),
               "test_arms": units.loc[units["split"] == "TEST", "arm"].tolist(), "lines": {}}
    for line in bm.LINES:
        regs, tgts, pos = regs_all[line], tgts_all[line], pos_all[line]
        regs["split"] = regs["arm"].map(asplit)
        regs["fold"] = regs["arm"].map(afold)
        regs["degree"] = regs["regulator"].map(pos[pos["pos_primary"]].groupby("regulator").size()).fillna(0)
        regs["degree_decile"] = bm.degree_deciles(regs, pos).values
        test_regs = set(regs.loc[regs["split"] == "TEST", "regulator"])
        lab = regs.loc[regs["split"] == "TEST", ["regulator", "own_class", "degree_decile"]]
        lab.to_csv(OUT / f"regulator_labels_{line}_TEST.tsv", sep="\t", index=False)
        out = regs.copy()
        out.loc[out["split"] == "TEST", ["own_class", "degree_decile"]] = np.nan
        cols = ["regulator", "pie_id", "chrom", "arm", "split", "fold", "own_class", "cn_tercile", "cn_stat",
                "degree_decile", "expr_context"]
        out["degree_decile"] = out["degree_decile"].astype("Int64")
        out["fold"] = out["fold"].astype("Int64").where(out["split"] == "TRAIN")
        out[cols].to_csv(OUT / f"universe_regulators_{line}.tsv", sep="\t", index=False)
        tgts[["target", "pie_id", "chrom", "arm", "mid", "bin", "expr_mean_s2", "expr_decile", "expr_hct116",
              "expr_hek293"]].to_csv(OUT / f"universe_targets_{line}.tsv", sep="\t", index=False)
        pos_te = pos[pos["regulator"].isin(test_regs)]
        pos_tr = pos[~pos["regulator"].isin(test_regs)]
        pos_tr.to_parquet(OUT / f"positives_{line}_TRAIN.parquet", index=False)
        pos_te.to_parquet(OUT / f"positives_{line}_TEST.parquet", index=False)
        ls = {}
        for sp in bm.SPLITS:
            r = regs[regs["split"] == sp]
            p = pos_tr if sp == "TRAIN" else pos_te
            ls[sp] = {"arms": int(r["arm"].nunique()), "regulators": len(r),
                      "pairs": int(sum(len(tgts) - (tgts["chrom"] == c).sum() for c in r["chrom"])),
                      **{f"n_{f}": int(p[f].sum()) for f in bm.POS_FLAGS}}
            if sp == "TRAIN":      # label-derived summaries are reported for TRAIN only (TEST stays sealed)
                ls[sp].update(frac_lfc_down=float((p.loc[p["pos_primary"], "lfc_sign"] < 0).mean()),
                              own_class=r["own_class"].value_counts().to_dict(),
                              regulators_with_ge1_positive=int(p["regulator"].nunique()),
                              median_on_target_lfc=float(r["own_lfc"].median()))
                manifest["lines"][line]["median_on_target_lfc_train"] = ls[sp]["median_on_target_lfc"]
                for f in (1, 2):
                    rf = r[r["fold"] == f]
                    ls[sp][f"fold{f}"] = {"arms": int(rf["arm"].nunique()), "regulators": len(rf),
                                          "positives": int(p.loc[p["regulator"].isin(rf["regulator"]),
                                                                 "pos_primary"].sum())}
        summary["lines"][line] = {"regulators": len(regs), "targets": len(tgts),
                                  "pairs": manifest["lines"][line]["n_pairs"], "splits": ls,
                                  "target_bins": int(tgts["bin"].nunique()),
                                  "excluded": {"regulators": manifest["lines"][line]["regulator_counts"],
                                               "targets": manifest["lines"][line]["target_counts"]},
                                  "n_cn_only_regulators": manifest["lines"][line]["n_cn_only_regulators"]}
    shared = set(regs_all["hct116"]["regulator"]) & set(regs_all["hek293t"]["regulator"])
    summary["shared_regulators"] = len(shared)
    manifest["outputs"] = {str(p.relative_to(ROOT)): bm.sha256_file(p) for p in sorted(OUT.glob("*"))
                           if p.name != "manifest.json"}
    manifest["outputs"].update({str(p.relative_to(ROOT)): bm.sha256_file(p) for p in sorted(SPLIT_DIR.glob("*.tsv"))})
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    (SPLIT_DIR / "benchmark_summary.json").write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
