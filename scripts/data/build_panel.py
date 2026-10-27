#!/usr/bin/env python
"""Build the frozen sample panel S2 and the frozen gene flags (DESIGN v3.1 §3; docs/INTERFACES.md,
"Sample panel and frozen flags"; data task T3).

S2 = DepMap 26Q3 cell lines (ModelType == "Cell Line"; no organoids) with matched default CN + RNA,
     excluding models with a non-empty Model.EngineeredModel, then one model per PatientID.
     Tie-break within a PatientID (explicit rule): (1) a model with 24Q4 PureCN values (sig24_Ploidy present),
     (2) the lowest ModelID (the earliest-registered, i.e. the parental line; DepMap has no other
     per-patient "default" flag).
in_S1 := the model has 24Q4 PureCN ploidy (S1 = PureCN-complete sensitivity panel; here flagged within S2).

Covariates are stored raw (encoding/pooling such as "levels >= 15" is the model's job):
  lineage, primary_disease, sex; ploidy, wgd (+ missing flags); msi_score (+ missing flag, msi_high =
  score > 20); TP53 (24Q4 rule: altered iff Damaging > 0 or Hotspot > 0; + missing flag); media_base (parsed
  from the FormulationID of the RNA ModelCondition, falling back to Model.FormulationID), media_formulation
  (raw), serum_free; growth_pattern; stranded, platform (RNA profile); PatientID; profile ids.

Frozen flags, one row per gene key (= DepMap column label, union of CN and expression columns):
  x = log2(clip(relCN, 2^-2, 2^3)) on S2 (DESIGN transform; |x| >= 1 is unaffected by the clip bounds).
  variable_cn: autosomal genes: >= 20 S2 lines with |x| >= 1. X genes (non-PAR): same test on sex-centred
    x (x minus the gene's median among S2 lines of the same Model.Sex; Unknown/missing sex form one group).
    Y, PAR and MT genes and genes without a placed chromosome: False (excluded by DESIGN §3).
  measured: has CN and expression column, all S2 values finite in both, SD > 0 in S2 for both.
Outputs: data/processed/panel_s2.tsv, data/processed/frozen_flags.tsv, data/processed/panel_summary.json
"""
from __future__ import annotations

import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data" / "processed"
MIN_LINES = 20
THRESH = 1.0
CLIP = (2.0 ** -2, 2.0 ** 3)
warnings.filterwarnings("ignore", category=RuntimeWarning)  # all-NaN (masked) CN genes


def media_base(f: str | float) -> str | None:
    if not isinstance(f, str) or not f.strip():
        return None
    s = f.strip()
    if re.match(r"^(DMEM\s*[:/]\s*F-?12|F-?12\s*[:/]\s*DMEM)", s, flags=re.I):
        return "DMEM/F12"
    if s.lower().startswith("tumor stem media"):
        return "TSM"
    return re.split(r"\s*\+\s*|\s*:\s*|\s*,\s*", s)[0].strip()


def build_panel() -> pd.DataFrame:
    m = pd.read_csv(PROC / "depmap_models.tsv", sep="\t", low_memory=False)
    tp = pd.read_csv(PROC / "depmap_tp53.tsv", sep="\t").set_index("ModelID")
    rt = pd.read_csv(PROC / "depmap_rna_technical.tsv", sep="\t").set_index("ModelID")
    c = m[m.matched & (m.model_ModelType == "Cell Line")].copy()
    n_matched_cl = len(c)
    c = c[c.model_EngineeredModel.isna()]
    n_after_eng = len(c)
    c["_pc"] = c.sig24_Ploidy.notna()
    c = c.sort_values(["model_PatientID", "_pc", "ModelID"], ascending=[True, False, True])
    assert c.model_PatientID.notna().all()
    s2 = c.drop_duplicates("model_PatientID").sort_values("ModelID")

    form = s2.mc_FormulationID.fillna(s2.model_FormulationID)
    sf = s2.mc_SerumFreeMedia.astype("boolean").fillna(s2.model_SerumFreeMedia.astype("boolean"))
    p = pd.DataFrame({
        "ModelID": s2.ModelID, "PatientID": s2.model_PatientID, "CellLineName": s2.model_CellLineName,
        "StrippedCellLineName": s2.model_StrippedCellLineName,
        "lineage": s2.model_OncotreeLineage, "primary_disease": s2.model_OncotreePrimaryDisease,
        "oncotree_code": s2.model_OncotreeCode, "sex": s2.model_Sex,
        "ploidy": s2.sig24_Ploidy, "ploidy_missing": s2.sig24_Ploidy.isna(),
        "wgd": s2.sig24_WGD, "wgd_missing": s2.sig24_WGD.isna(),
        "msi_score": s2.sig24_MSIScore, "msi_missing": s2.sig24_MSIScore.isna(),
        "msi_high": (s2.sig24_MSIScore > 20).where(s2.sig24_MSIScore.notna()),
        "sig24_source": s2.sig24_source,
    })
    t = tp.reindex(s2.ModelID)
    p["tp53_status"] = t.tp53_status.to_numpy()
    p["tp53_altered"] = (t.tp53_status == "mutant").where(t.tp53_status != "unknown").to_numpy()
    p["tp53_missing"] = (t.tp53_status == "unknown").to_numpy()
    p["media_base"] = form.map(media_base).to_numpy()
    p["media_formulation"] = form.to_numpy()
    p["serum_free"] = sf.to_numpy()
    p["growth_pattern"] = s2.model_GrowthPattern.to_numpy()
    r = rt.reindex(s2.ModelID)
    p["stranded"] = r.rna_Stranded.to_numpy()
    p["platform"] = r.rna_SequencingPlatform.to_numpy()
    p["in_S1"] = s2.sig24_Ploidy.notna().to_numpy()
    p["cn_seq_id"] = s2.cn_seq_id.to_numpy()
    p["rna_seq_id"] = s2.rna_seq_id.to_numpy()
    p["same_mc"] = s2.same_mc.to_numpy()
    p.attrs["counts"] = {"matched_cell_lines": n_matched_cl, "after_engineered_exclusion": n_after_eng,
                         "S2": len(p), "dropped_patient_duplicates": n_after_eng - len(p)}
    return p.reset_index(drop=True)


def build_flags(panel: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    ann = pd.read_csv(PROC / "gene_annotation.tsv", sep="\t", dtype={"chrom": str}).set_index("depmap_column")
    cn = pl.read_parquet(PROC / "depmap_cn_matched.parquet")
    ex = pl.read_parquet(PROC / "depmap_expr_matched.parquet")
    keep = panel.ModelID.tolist()
    cn = cn.filter(pl.col("ModelID").is_in(keep)).sort("ModelID")
    ex = ex.filter(pl.col("ModelID").is_in(keep)).sort("ModelID")
    assert cn["ModelID"].to_list() == ex["ModelID"].to_list() == sorted(keep)
    sex = panel.set_index("ModelID").sex.reindex(cn["ModelID"].to_list()).fillna("Unknown").to_numpy()
    cn_cols = [c for c in cn.columns if c != "ModelID"]
    ex_cols = [c for c in ex.columns if c != "ModelID"]
    X = cn.select(cn_cols).to_numpy().astype(np.float64)
    E = ex.select(ex_cols).to_numpy().astype(np.float64)
    L = np.log2(np.clip(X, *CLIP))
    a = ann.reindex(cn_cols)
    is_x = ((a.chrom == "X") & ~a.is_par.astype(bool)).to_numpy()
    Lc = L.copy()
    for s in np.unique(sex):
        rows = sex == s
        Lc[np.ix_(rows, is_x)] = L[np.ix_(rows, is_x)] - np.nanmedian(L[np.ix_(rows, is_x)], axis=0)
    n_var = np.sum(np.abs(np.where(is_x[None, :], Lc, L)) >= THRESH, axis=0)
    n_var_raw_x = np.sum(np.abs(L) >= THRESH, axis=0)
    # robustness: X test excluding Unknown-sex lines
    known = np.isin(sex, ["Male", "Female"])
    n_var_x_known = np.sum(np.abs(Lc[known]) >= THRESH, axis=0)
    cn_ok = np.isfinite(X).all(0) & (np.nanstd(X, 0) > 0)
    ex_ok = np.isfinite(E).all(0) & (np.nanstd(E, 0) > 0)

    keys = sorted(set(cn_cols) | set(ex_cols))
    f = pd.DataFrame({"gene_key": keys})
    an = ann.reindex(keys)
    f["symbol"] = an.symbol.to_numpy()
    f["chrom"] = an.chrom.to_numpy()
    f["arm"] = an.arm.to_numpy()
    f["chrom_class"] = an.chrom_class.to_numpy()
    f["has_cn"] = f.gene_key.isin(cn_cols)
    f["has_expr"] = f.gene_key.isin(ex_cols)
    f["is_par"] = an.is_par.fillna(False).astype(bool).to_numpy()
    f["is_y"] = (an.chrom == "Y").to_numpy()
    f["is_mt"] = (an.chrom == "MT").to_numpy()
    ci = {c: i for i, c in enumerate(cn_cols)}
    ei = {c: i for i, c in enumerate(ex_cols)}
    f["n_lines_abs_log2cn_ge1"] = [int(n_var[ci[k]]) if k in ci else np.nan for k in keys]
    eligible = f.chrom.isin([str(i) for i in range(1, 23)] + ["X"]) & ~f.is_par & f.has_cn
    f["cn_test"] = np.where(~eligible, "excluded", np.where(f.chrom == "X", "sex_centred", "raw"))
    f["variable_cn"] = eligible & (f.n_lines_abs_log2cn_ge1 >= MIN_LINES)
    f["measured"] = [bool(k in ci and k in ei and cn_ok[ci[k]] and ex_ok[ei[k]]) for k in keys]

    xs = f[f.cn_test == "sex_centred"]
    s = {
        "n_gene_keys": len(f), "n_samples": len(keep),
        "variable_cn": {"autosomal": int((f.variable_cn & (f.cn_test == "raw")).sum()),
                        "X_sex_centred": int((f.variable_cn & (f.cn_test == "sex_centred")).sum()),
                        "total": int(f.variable_cn.sum())},
        "X_eligible": int(len(xs)),
        "X_variable_if_not_sex_centred": int(sum(n_var_raw_x[ci[k]] >= MIN_LINES for k in xs.gene_key)),
        "X_variable_excluding_unknown_sex": int(sum(n_var_x_known[ci[k]] >= MIN_LINES for k in xs.gene_key)),
        "autosomal_eligible": int((f.cn_test == "raw").sum()),
        "measured": int(f.measured.sum()),
        "measured_excl_y_par_mt": int((f.measured & ~f.is_y & ~f.is_par & ~f.is_mt).sum()),
        "measured_by_chrom_class": f[f.measured].chrom_class.value_counts().to_dict(),
        "has_cn_and_expr": int((f.has_cn & f.has_expr).sum()),
        "not_measured_but_both": int((f.has_cn & f.has_expr & ~f.measured).sum()),
        "variable_and_measured": int((f.variable_cn & f.measured).sum()),
        "n_lines_quantiles": {t: f[(f.cn_test == t) & (f.n_lines_abs_log2cn_ge1 > 0)].n_lines_abs_log2cn_ge1
                              .quantile([0.05, 0.25, 0.5, 0.75]).to_dict() for t in ["raw", "sex_centred"]},
        "female_lines_with_median_X_log2cn_lt_-0.7": float(np.mean(
            np.nanmedian(L[np.ix_(sex == "Female", is_x)], axis=1) < -0.7)),
        "sex_groups": {str(k): int(v) for k, v in zip(*np.unique(sex, return_counts=True))},
    }
    return f, s


def main() -> int:
    panel = build_panel()
    counts = panel.attrs["counts"]
    panel.to_csv(PROC / "panel_s2.tsv", sep="\t", index=False)
    flags, fs = build_flags(panel)
    flags.to_csv(PROC / "frozen_flags.tsv", sep="\t", index=False)
    summ = {"panel": counts | {
        "in_S1": int(panel.in_S1.sum()), "lineages": int(panel.lineage.nunique()),
        "ploidy_missing": int(panel.ploidy_missing.sum()), "msi_missing": int(panel.msi_missing.sum()),
        "msi_high": int((panel.msi_high == True).sum()),  # noqa: E712
        "tp53": panel.tp53_status.value_counts().to_dict(),
        "sex": panel.sex.fillna("missing").value_counts().to_dict(),
        "media_base_levels": int(panel.media_base.nunique()), "media_base_missing": int(panel.media_base.isna().sum()),
        "media_base_levels_ge15": {k: int(v) for k, v in panel.media_base.value_counts().items() if v >= 15},
        "serum_free": panel.serum_free.astype(str).value_counts().to_dict(),
        "growth_pattern": panel.growth_pattern.value_counts().to_dict(),
        "stranded": panel.stranded.astype(str).value_counts().to_dict(),
        "platform_missing": int(panel.platform.isna().sum()),
        "benchmark_lines": {n: bool((panel.StrippedCellLineName == n).any()) for n in ["HCT116", "HEK293"]},
    }, "flags": fs}
    (PROC / "panel_summary.json").write_text(json.dumps(summ, indent=1, default=str))
    print(json.dumps(summ, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
