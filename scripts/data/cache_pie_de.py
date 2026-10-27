#!/usr/bin/env python
"""Build compact caches of the PIE X-Atlas/Orion DE tables so later steps never rescan the 5 GB parquets.

Input : data/external/pie_x_atlas_orion/de/{hct116,hek293t}.parquet (HF commit abfa3da0...)
Output (data/processed/pie/):
  vocab_targets.tsv, vocab_features.tsv   global integer codes shared by both cell lines
  tested_<cell>.parquet   ALL rows (one per tested (target, feature) pair), compact:
                          t (uint16 target code), f (uint16 feature code), lfc (float32 log2FC),
                          mlog10_p, mlog10_padj (float32 -log10 of p_value / p_adj; p == 0 -> +inf)
                          "tested" := row present in the DE table (see docs/DATA_PROFILE.md).
  loose_<cell>.parquet    rows with p_value < LOOSE_P OR p_adj < LOOSE_PADJ at full precision with string ids
                          (target, feature, log2_fold_change, p_value, p_adj, Ueffect)
  targets_<cell>.tsv      per target: target_ncells, ref_ncells, n_features tested
  features_<cell>.tsv     per feature: ref_mean (control mean expression), n_targets tested
Also checks: (target, feature) uniqueness, per-target constancy of target_ncells, per-feature
constancy of ref_mean, file size vs expected.

Usage: .venv/bin/python scripts/data/cache_pie_de.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data" / "external" / "pie_x_atlas_orion" / "de"
OUT = ROOT / "data" / "processed" / "pie"
CELLS = {"hct116": 3395023058, "hek293t": 1820568723}  # expected byte sizes (HF commit abfa3da0)
LOOSE_P = 1e-3
LOOSE_PADJ = 0.1  # guarantees every p_adj < 0.1 row is present at full precision
BATCH = 4_000_000


def mlog10(x: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore"):
        return (-np.log10(x)).astype(np.float32)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    for cell, size in CELLS.items():
        p = SRC / f"{cell}.parquet"
        if p.stat().st_size != size:
            raise SystemExit(f"{p} has size {p.stat().st_size}, expected {size} (download incomplete?)")

    # pass 1: global vocabularies
    targets, features = set(), set()
    for cell in CELLS:
        lf = pl.scan_parquet(SRC / f"{cell}.parquet")
        targets |= set(lf.select(pl.col("target").unique()).collect(engine="streaming")["target"].to_list())
        features |= set(lf.select(pl.col("feature").unique()).collect(engine="streaming")["feature"].to_list())
    tv = sorted(targets)
    fv = sorted(features)
    assert len(tv) < 65535 and len(fv) < 65535
    pd.DataFrame({"code": range(len(tv)), "target": tv}).to_csv(OUT / "vocab_targets.tsv", sep="\t", index=False)
    pd.DataFrame({"code": range(len(fv)), "feature": fv}).to_csv(OUT / "vocab_features.tsv", sep="\t", index=False)
    print(f"vocab: {len(tv)} targets, {len(fv)} features", file=sys.stderr)
    tcat = pd.CategoricalDtype(tv)
    fcat = pd.CategoricalDtype(fv)

    # pass 2: compact tables
    for cell in CELLS:
        pf = pq.ParquetFile(SRC / f"{cell}.parquet")
        cols = ["target", "feature", "log2_fold_change", "p_value", "p_adj", "Ueffect",
                "target_ncells", "ref_ncells", "ref_mean"]
        tested_writer = None
        loose_parts = []
        tinfo = []
        finfo = []
        n = 0
        for batch in pf.iter_batches(batch_size=BATCH, columns=cols):
            df = batch.to_pandas()
            n += len(df)
            t = df["target"].astype(tcat).cat.codes.to_numpy().astype(np.uint16)
            f = df["feature"].astype(fcat).cat.codes.to_numpy().astype(np.uint16)
            tbl = pa.table({
                "t": t, "f": f,
                "lfc": df["log2_fold_change"].to_numpy(np.float32),
                "mlog10_p": mlog10(df["p_value"].to_numpy()),
                "mlog10_padj": mlog10(df["p_adj"].to_numpy()),
            })
            if tested_writer is None:
                tested_writer = pq.ParquetWriter(OUT / f"tested_{cell}.parquet", tbl.schema, compression="zstd")
            tested_writer.write_table(tbl)
            loose_parts.append(df.loc[(df["p_value"] < LOOSE_P) | (df["p_adj"] < LOOSE_PADJ),
                                      ["target", "feature", "log2_fold_change", "p_value", "p_adj", "Ueffect"]])
            tinfo.append(df.groupby("target", sort=False).agg(
                target_ncells_min=("target_ncells", "min"), target_ncells_max=("target_ncells", "max"),
                ref_ncells=("ref_ncells", "max"), n_features=("feature", "size")))
            finfo.append(df.groupby("feature", sort=False).agg(
                ref_mean_min=("ref_mean", "min"), ref_mean_max=("ref_mean", "max"), n_targets=("target", "size")))
        tested_writer.close()
        if n != pf.metadata.num_rows:
            raise SystemExit("row count mismatch")

        ti = pd.concat(tinfo).groupby(level=0).agg(
            target_ncells_min=("target_ncells_min", "min"), target_ncells_max=("target_ncells_max", "max"),
            ref_ncells=("ref_ncells", "max"), n_features=("n_features", "sum"))
        fi = pd.concat(finfo).groupby(level=0).agg(
            ref_mean_min=("ref_mean_min", "min"), ref_mean_max=("ref_mean_max", "max"), n_targets=("n_targets", "sum"))
        if (ti.target_ncells_min != ti.target_ncells_max).any():
            print(f"WARNING {cell}: target_ncells not constant within target", file=sys.stderr)
        rel = (fi.ref_mean_max - fi.ref_mean_min) / fi.ref_mean_max.where(fi.ref_mean_max > 0, 1)
        if (rel > 1e-6).any():
            print(f"WARNING {cell}: ref_mean varies within feature for {(rel > 1e-6).sum()} features", file=sys.stderr)
        ti = ti.rename(columns={"target_ncells_max": "target_ncells"}).drop(columns="target_ncells_min")
        fi = fi.rename(columns={"ref_mean_max": "ref_mean"}).drop(columns="ref_mean_min")
        ti.index.name = "target"
        fi.index.name = "feature"
        ti.to_csv(OUT / f"targets_{cell}.tsv", sep="\t")
        fi.to_csv(OUT / f"features_{cell}.tsv", sep="\t")

        loose = pd.concat(loose_parts, ignore_index=True)
        loose.to_parquet(OUT / f"loose_{cell}.parquet", compression="zstd", index=False)

        # uniqueness of (target, feature)
        dup = (pl.scan_parquet(OUT / f"tested_{cell}.parquet").group_by(["t", "f"]).len()
               .filter(pl.col("len") > 1).collect(engine="streaming").height)
        print(f"{cell}: {n} rows, {len(ti)} targets, {len(fi)} features, duplicate pairs {dup}, "
              f"loose(p<{LOOSE_P} or padj<{LOOSE_PADJ}) rows {len(loose)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
