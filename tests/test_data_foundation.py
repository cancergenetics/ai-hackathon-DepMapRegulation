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


# ---------------- DESIGN v2 data tasks (T1, T2/E4, T4) ----------------
def test_par_flags_and_x_copy_link():
    a = pd.read_csv(PROC / "gene_annotation.tsv", sep="\t", dtype={"chrom": str})
    py = a[a.depmap_column.str.endswith("_PAR_Y")]
    assert len(py) == 3 and py.is_par.all() and py.par_y_copy.all() and (py.chrom == "Y").all()
    assert a[a.is_par].par_region.notna().all()
    p = pd.read_csv(PROC / "pie_gene_map.tsv", sep="\t", keep_default_na=False, na_values=[""])
    cn = p.depmap_cn_column.dropna()
    assert not cn.str.contains(";").any() and not cn.str.endswith("_PAR_Y").any()
    assert p.set_index("pie_id").loc["IL9R", "depmap_cn_column"] == "IL9R (3581)"


def test_tp53_table():
    t = pd.read_csv(PROC / "depmap_tp53.tsv", sep="\t")
    m = pd.read_csv(PROC / "depmap_models.tsv", sep="\t", low_memory=False)
    assert t.ModelID.is_unique and set(t.ModelID) == set(m.ModelID)
    assert set(t.tp53_status) <= {"mutant", "wt", "unknown"}
    g = t[t.genotyped_24Q4]
    assert ((g.tp53_status == "mutant") == ((g.tp53_damaging > 0) | (g.tp53_hotspot > 0))).all()
    assert (t[~t.genotyped_24Q4].tp53_status == "unknown").all()
    st = t.merge(m[["ModelID", "model_StrippedCellLineName"]]).set_index("model_StrippedCellLineName").tp53_status
    for line in ["HCT116", "MCF7", "A549", "U2OS"]:          # TP53 wild-type controls
        assert st[line] == "wt", line
    for line in ["SKBR3", "SW480", "HT29", "BT474", "T47D", "MDAMB231"]:  # TP53-mutant controls
        assert st[line] == "mutant", line


def test_rna_technical_table():
    r = pd.read_csv(PROC / "depmap_rna_technical.tsv", sep="\t")
    m = pd.read_csv(PROC / "depmap_models.tsv", sep="\t", low_memory=False)
    assert r.ModelID.is_unique and set(r.ModelID) == set(m[m.has_rna].ModelID)
    assert r.rna_Stranded.notna().all()
    assert (r.set_index("ModelID").rna_seq_id == m.set_index("ModelID").rna_seq_id.reindex(r.ModelID)).all()


def test_patient_id_present_for_cell_lines():
    m = pd.read_csv(PROC / "depmap_models.tsv", sep="\t", low_memory=False)
    assert m[m.matched & (m.model_ModelType == "Cell Line")].model_PatientID.notna().all()


# ---------------- T3: frozen panel S2 and gene flags ----------------
def test_panel_s2():
    p = pd.read_csv(PROC / "panel_s2.tsv", sep="\t")
    m = pd.read_csv(PROC / "depmap_models.tsv", sep="\t", low_memory=False).set_index("ModelID")
    assert len(p) == 1436  # DESIGN v3.1 §3 (round-2 data check)
    assert p.ModelID.is_unique and p.PatientID.is_unique
    mm = m.loc[p.ModelID]
    assert mm.matched.all() and (mm.model_ModelType == "Cell Line").all() and mm.model_EngineeredModel.isna().all()
    assert "HCT116" in set(p.StrippedCellLineName)
    assert (p.in_S1 == p.ploidy.notna()).all() and (p.ploidy_missing == p.ploidy.isna()).all()
    assert (p.tp53_missing == (p.tp53_status == "unknown")).all() and int((~p.tp53_missing).sum()) == 1255
    # dedup preference: if any eligible model of a patient has PureCN, the chosen one does
    elig = m[m.matched & (m.model_ModelType == "Cell Line") & m.model_EngineeredModel.isna()]
    has_pc = elig.groupby("model_PatientID").sig24_Ploidy.apply(lambda s: s.notna().any())
    chosen = p.set_index("PatientID").in_S1
    assert (chosen[has_pc[has_pc].index]).all()


def test_frozen_flags():
    f = pd.read_csv(PROC / "frozen_flags.tsv", sep="\t", dtype={"chrom": str}).set_index("gene_key")
    a = pd.read_csv(PROC / "gene_annotation.tsv", sep="\t")
    assert f.index.is_unique and set(f.index) == set(a.depmap_column)
    assert not f[f.is_y | f.is_par | f.is_mt].variable_cn.any()
    assert (f[f.measured].has_cn & f[f.measured].has_expr).all()
    for g in ["ERBB2 (2064)", "MYC (4609)", "CDKN2A (1029)", "EGFR (1956)"]:  # known recurrent CN events
        assert f.loc[g, "variable_cn"] and f.loc[g, "measured"], g
    # recompute one autosomal and one X gene count from the matrices
    p = pd.read_csv(PROC / "panel_s2.tsv", sep="\t")
    xg = f[(f.cn_test == "sex_centred") & (f.n_lines_abs_log2cn_ge1 > 0)].index[0]
    cn = pl.read_parquet(PROC / "depmap_cn_matched.parquet", columns=["ModelID", "MYC (4609)", xg]).to_pandas()
    cn = cn[cn.ModelID.isin(p.ModelID)].set_index("ModelID")
    L = np.log2(np.clip(cn.astype(float), 0.25, 8))
    assert int((L["MYC (4609)"].abs() >= 1).sum()) == f.loc["MYC (4609)", "n_lines_abs_log2cn_ge1"]
    sex = p.set_index("ModelID").sex.reindex(L.index).fillna("Unknown")
    xc = L[xg] - L[xg].groupby(sex).transform("median")
    assert int((xc.abs() >= 1).sum()) == f.loc[xg, "n_lines_abs_log2cn_ge1"]
