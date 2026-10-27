"""Sanity checks on the data foundation (scripts/data/*). Run after scripts/data/run_all.sh:
    .venv/bin/python -m pytest -q tests/test_data_foundation.py
Positive controls: known amplicons / deletions, chromosome assignment of landmark genes, on-target
CRISPRi knockdown. Structural checks: one row per model, aligned matrices, one-to-one gene maps,
BH applied within target. Skipped if outputs are absent (data/ is not committed).
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pyarrow.parquet as pq
import pytest

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
pytestmark = pytest.mark.skipif(not (PROC / "gene_annotation.tsv").exists(), reason="processed data absent")
sys.path.insert(0, str(ROOT / "scripts" / "data"))


def test_bsd_sum_toy(tmp_path):
    from download_annotation import bsd_sum
    p = tmp_path / "x.bin"
    p.write_bytes(b"hello world\n" * 1000)
    try:  # macOS/BSD `sum` or GNU `sum -r` implement the same algorithm
        out = subprocess.run(["sum", p], capture_output=True, text=True, check=True).stdout.split()
    except Exception:
        pytest.skip("no sum binary")
    assert bsd_sum(p) == (int(out[0]), int(out[1]))


def test_landmark_gene_locations():
    a = pd.read_csv(PROC / "gene_annotation.tsv", sep="\t", dtype={"chrom": str}).set_index("symbol")
    exp = {"ERBB2": ("17", "q"), "MYC": ("8", "q"), "TP53": ("17", "p"), "CDKN2A": ("9", "p"), "XIST": ("X", "q"),
           "SRY": ("Y", "p"), "EGFR": ("7", "p")}
    for g, (c, arm) in exp.items():
        if g in a.index:
            r = a.loc[g]
            r = r.iloc[0] if isinstance(r, pd.DataFrame) else r
            assert (r.chrom, r.arm) == (c, arm), g
    assert a.loc["MT-CO1"].chrom == "MT"


def test_annotation_structure():
    a = pd.read_csv(PROC / "gene_annotation.tsv", sep="\t", dtype={"chrom": str})
    assert a.depmap_column.is_unique
    assert not a.hgnc_dup_in_depmap.any()
    # every expression gene is placed on a primary chromosome
    assert (a[a.in_expr].chrom_class != "unplaced").all()


def test_pie_map_one_to_one():
    p = pd.read_csv(PROC / "pie_gene_map.tsv", sep="\t", keep_default_na=False, na_values=[""])
    k = p[p.depmap_expr_column.notna()]
    assert k.depmap_expr_column.is_unique
    k = p[p.depmap_cn_column.notna()]
    assert k.depmap_cn_column.is_unique


def test_models_and_matrices_aligned():
    m = pd.read_csv(PROC / "depmap_models.tsv", sep="\t", low_memory=False)
    assert m.ModelID.is_unique
    mm = m[m.matched].ModelID.tolist()
    for f in ["depmap_cn_matched.parquet", "depmap_expr_matched.parquet"]:
        assert pl.read_parquet(PROC / f, columns=["ModelID"])["ModelID"].to_list() == mm


def test_known_copy_number_events():
    m = pd.read_csv(PROC / "depmap_models.tsv", sep="\t", low_memory=False)
    m = m[m.matched].set_index("model_StrippedCellLineName")
    cn = pl.read_parquet(PROC / "depmap_cn_matched.parquet",
                         columns=["ModelID", "ERBB2 (2064)", "CDKN2A (1029)", "MYC (4609)"]).to_pandas().set_index("ModelID")
    def v(line, g):
        return float(cn.loc[m.loc[line, "ModelID"], g])
    assert v("SKBR3", "ERBB2 (2064)") > 5 and v("BT474", "ERBB2 (2064)") > 5   # ERBB2 amplified
    assert v("MCF7", "ERBB2 (2064)") < 1.5                                   # ERBB2 not amplified
    assert v("COLO320", "MYC (4609)") > 20                                   # MYC ecDNA
    assert v("A549", "CDKN2A (1029)") < 0.1 and v("U87MG", "CDKN2A (1029)") < 0.1  # homozygous deletion
    # relative CN is linear and centred near 1
    prof = json.loads((PROC / "depmap_profile.json").read_text())
    assert 0.95 < prof["cn_values"]["quantiles"]["0.5"] < 1.1
    assert prof["cn_vs_absolute_24Q4"]["corr_ploidy_vs_median_abs_over_rel"] > 0.9


def test_benchmark_on_target_and_bh():
    b = json.loads((PROC / "benchmark_profile.json").read_text())
    for cell, r in b["per_cell"].items():
        assert r["on_target"]["frac_lfc_lt_0"] > 0.95, cell
        n = pq.ParquetFile(ROOT / "data/external/pie_x_atlas_orion/de" / f"{cell}.parquet").metadata.num_rows
        assert pq.ParquetFile(PROC / "pie" / f"tested_{cell}.parquet").metadata.num_rows == n
        for tg, chk in b["bh_scope_check"][cell].items():
            assert chk["max_abs_diff_within_target_BH"] < 1e-9, (cell, tg)


def test_loose_subset_complete_for_padj_0_1():
    for cell in ["hct116", "hek293t"]:
        n_loose = pl.scan_parquet(PROC / "pie" / f"loose_{cell}.parquet").filter(pl.col("p_adj") < 0.05).select(pl.len()).collect().item()
        n_all = pl.scan_parquet(PROC / "pie" / f"tested_{cell}.parquet").filter(
            pl.col("mlog10_padj") > np.float32(-np.log10(0.05))).select(pl.len()).collect().item()
        assert abs(n_loose - n_all) <= max(5, 1e-4 * n_all), (cell, n_loose, n_all)  # float32 rounding at threshold
