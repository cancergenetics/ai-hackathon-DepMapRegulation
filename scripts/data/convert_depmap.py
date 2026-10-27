#!/usr/bin/env python
"""Convert the DepMap 26Q3 gene-level CN and expression CSVs (profile-level rows) to parquet.

The CSVs have one row per sequencing profile with metadata columns
(SequencingID, ModelConditionID, ModelID, IsDefaultEntryForMC, IsDefaultEntryForModel) followed by
gene columns "SYMBOL (ENTREZ)". The first (unnamed) column is a 0..N-1 pandas row number and is
dropped after checking that. Values are stored as float32 (CSV precision is ~6-7 significant
digits, so no information is lost in practice; the check below quantifies max abs rounding error).

Usage: .venv/bin/python scripts/data/convert_depmap.py
Writes data/interim/depmap_cn_wgs.parquet, data/interim/depmap_expr_tpm_logp1.parquet
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "starting_data"
OUT = ROOT / "data" / "interim"
META = ["SequencingID", "ModelConditionID", "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel"]
FILES = {
    "OmicsCNGeneWGS.csv": "depmap_cn_wgs.parquet",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv": "depmap_expr_tpm_logp1.parquet",
}


def convert(src: Path, dest: Path) -> None:
    header = pl.read_csv(src, n_rows=0).columns
    idx_col = header[0]
    genes = [c for c in header if c not in META and c != idx_col]
    schema = {c: pl.Utf8 for c in [idx_col, *META]} | {g: pl.Float64 for g in genes}
    df = pl.read_csv(src, schema_overrides=schema, infer_schema=False)
    if not (df[idx_col].cast(pl.Int64) == pl.int_range(df.height, eager=True)).all():
        raise SystemExit(f"{src.name}: unexpected leading index column")
    x64 = df.select(genes).to_numpy()
    x32 = x64.astype(np.float32)
    err = np.nanmax(np.abs(x64 - x32.astype(np.float64)))
    out = pl.concat([df.select(META), pl.DataFrame(x32, schema=genes)], how="horizontal")
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.write_parquet(dest, compression="zstd", statistics=False)
    print(f"{src.name}: {out.height} rows x {len(genes)} genes -> {dest.name}; float32 max abs err {err:.2e}",
          file=sys.stderr)


def main() -> int:
    for s, d in FILES.items():
        convert(SRC / s, OUT / d)
    return 0


if __name__ == "__main__":
    sys.exit(main())
