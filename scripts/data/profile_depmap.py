#!/usr/bin/env python
"""Profile the DepMap 26Q3 CN/expression foundation. Writes data/processed/depmap_profile.json
(all numbers quoted in docs/DATA_PROFILE.md) and prints a readable summary.

Requires: convert_depmap.py, build_gene_annotation.py, build_depmap_tables.py outputs.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
SD = ROOT / "starting_data"
EXT = ROOT / "data" / "external"
PROC = ROOT / "data" / "processed"
INTERIM = ROOT / "data" / "interim"
META = ["SequencingID", "ModelConditionID", "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel"]
CN_FLOOR = 2 ** -5  # same floor as build_depmap_tables.py
warnings.filterwarnings("ignore", category=RuntimeWarning)


def q(x, qs=(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)) -> dict:
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    return {str(k): float(v) for k, v in zip(qs, np.quantile(x, qs))} if len(x) else {}


def vc(s: pd.Series, n: int = 15) -> dict:
    return {str(k): int(v) for k, v in s.value_counts(dropna=False).head(n).items()}


def main() -> int:
    R: dict = {}
    models = pd.read_csv(PROC / "depmap_models.tsv", sep="\t", low_memory=False)
    m = models[models.matched].reset_index(drop=True)
    ann = pd.read_csv(PROC / "gene_annotation.tsv", sep="\t", dtype={"chrom": str, "entrez_id": str})
    pie = pd.read_csv(PROC / "pie_gene_map.tsv", sep="\t", keep_default_na=False, na_values=[""],
                      dtype={"chrom": str})

    # ---------------- samples ----------------
    cn_meta = pl.read_parquet(INTERIM / "depmap_cn_wgs.parquet", columns=META).to_pandas()
    ex_meta = pl.read_parquet(INTERIM / "depmap_expr_tpm_logp1.parquet", columns=META).to_pandas()
    prof = pd.read_csv(SD / "OmicsProfiles.csv")
    R["samples"] = {
        "cn_rows": len(cn_meta), "cn_models": int(cn_meta.ModelID.nunique()),
        "cn_default_rows": int((cn_meta.IsDefaultEntryForModel == "Yes").sum()),
        "cn_models_with_multiple_rows": int((cn_meta.ModelID.value_counts() > 1).sum()),
        "rna_rows": len(ex_meta), "rna_models": int(ex_meta.ModelID.nunique()),
        "rna_default_rows": int((ex_meta.IsDefaultEntryForModel == "Yes").sum()),
        "rna_models_with_multiple_rows": int((ex_meta.ModelID.value_counts() > 1).sum()),
        "matched_models": int(models.matched.sum()),
        "cn_only_models": int((models.has_cn & ~models.has_rna).sum()),
        "rna_only_models": int((models.has_rna & ~models.has_cn).sum()),
        "matched_same_model_condition": int((m.same_mc == True).sum()),  # noqa: E712
        "matched_model_type": vc(m.model_ModelType),
        "matched_depmap_model_type_n": int(m.model_DepmapModelType.nunique()),
        "matched_any_rows_models": len(set(cn_meta.ModelID) & set(ex_meta.ModelID)),
        "omics_profiles_datatype": vc(prof.DataType),
        "rna_stranded_matched": vc(m.rna_stranded),
        "cn_seq_year_matched": vc(pd.to_datetime(m.cn_seq_date, errors="coerce").dt.year, 20),
        "matched_model_sex": vc(m.model_Sex),
        "matched_engineered": int(m.model_EngineeredModel.notna().sum()),
    }
    names = models.set_index("model_StrippedCellLineName")
    R["samples"]["benchmark_lines_in_depmap"] = {
        nm: (dict(ModelID=names.loc[nm, "ModelID"], matched=bool(names.loc[nm, "matched"]))
             if nm in names.index else None) for nm in ["HCT116", "HEK293T", "HEKTE", "HEK293"]}

    # ---------------- genes ----------------
    cn_cols = [c for c in pl.read_parquet_schema(PROC / "depmap_cn_matched.parquet") if c != "ModelID"]
    ex_cols = [c for c in pl.read_parquet_schema(PROC / "depmap_expr_matched.parquet") if c != "ModelID"]
    both = sorted(set(cn_cols) & set(ex_cols))
    a = ann.set_index("depmap_column")
    cn_only = a.loc[sorted(set(cn_cols) - set(ex_cols))]
    ex_only = a.loc[sorted(set(ex_cols) - set(cn_cols))]
    R["genes"] = {
        "cn_columns": len(cn_cols), "expr_columns": len(ex_cols), "both_columns": len(both),
        "both_by_hgnc": len(set(a.loc[cn_cols].hgnc_id.dropna()) & set(a.loc[ex_cols].hgnc_id.dropna())),
        "cn_only": len(cn_only), "expr_only": len(ex_only),
        "cn_only_locus_group": vc(cn_only.locus_group.fillna("no HGNC (bare Ensembl id)")),
        "cn_only_readthrough": int(cn_only.locus_type.fillna("").str.contains("readthrough").sum()),
        "expr_only_locus_type": vc(ex_only.locus_type.fillna("none")),
        "both_chrom_class": vc(a.loc[both].chrom_class),
    }

    # ---------------- CN values ----------------
    cn = pl.read_parquet(PROC / "depmap_cn_matched.parquet")
    ex = pl.read_parquet(PROC / "depmap_expr_matched.parquet")
    assert cn["ModelID"].to_list() == ex["ModelID"].to_list() == m.ModelID.tolist()
    X = cn.select(cn_cols).to_numpy().astype(np.float64)
    E = ex.select(ex_cols).to_numpy().astype(np.float64)
    v = X[~np.isnan(X)]
    L = np.log2(np.maximum(X, CN_FLOOR))
    R["cn_values"] = {
        "quantiles": q(v, (0, 0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999, 1)),
        "frac_negative": float((v < 0).mean()), "frac_zero": float((v == 0).mean()),
        "frac_lt_0.25": float((v < 0.25).mean()), "frac_gt_2": float((v > 2).mean()),
        "frac_gt_4": float((v > 4).mean()),
        "per_model_median": q(np.nanmedian(X, 1)), "per_model_mean": q(np.nanmean(X, 1)),
        "frac_below_floor": float((v < CN_FLOOR).mean()),
    }
    ev = E[~np.isnan(E)]
    R["expr_values"] = {"quantiles": q(ev), "frac_zero": float((ev == 0).mean()),
                        "per_model_frac_zero": q((E == 0).mean(1))}

    # relation to 24Q4 PureCN absolute CN / ploidy (models whose 24Q4 signature is the same profile)
    absg = pd.read_csv(EXT / "depmap_24Q4" / "OmicsAbsoluteCNGene.csv", index_col=0)
    cnpd = pd.DataFrame(X, index=m.ModelID, columns=cn_cols)
    common_g = sorted(set(absg.columns) & set(cn_cols))
    sp = m.set_index("ModelID")
    sel = [mid for mid in m.ModelID if mid in absg.index and sp.loc[mid, "sig24_source"] == "same_profile"]
    rows = []
    for mid in sel:
        r_ = cnpd.loc[mid, common_g].to_numpy(float)
        a_ = absg.loc[mid, common_g].to_numpy(float)
        ok = ~np.isnan(r_) & ~np.isnan(a_) & (r_ > 0.05)
        if ok.sum() < 1000:
            continue
        rows.append(dict(ModelID=mid, ploidy=sp.loc[mid, "sig24_Ploidy"],
                         median_abs_over_rel=float(np.median(a_[ok] / r_[ok])),
                         spearman_rel_abs=float(stats.spearmanr(r_[ok], a_[ok]).statistic),
                         pearson_rel_abs_over_ploidy=float(np.corrcoef(r_[ok], a_[ok] / sp.loc[mid, "sig24_Ploidy"])[0, 1]),
                         pearson_log2rel_vs_log2abs=float(np.corrcoef(np.log2(r_[ok]), np.log2(np.maximum(a_[ok], 0.5)))[0, 1]),
                         mean_rel=float(np.nanmean(r_))))
    rel = pd.DataFrame(rows)
    R["cn_vs_absolute_24Q4"] = {
        "n_models": len(rel), "n_common_genes": len(common_g),
        "median_spearman_rel_abs": float(rel.spearman_rel_abs.median()),
        "median_pearson_rel_vs_abs_over_ploidy": float(rel.pearson_rel_abs_over_ploidy.median()),
        "median_pearson_log2rel_vs_log2abs": float(rel.pearson_log2rel_vs_log2abs.median()),
        "corr_ploidy_vs_median_abs_over_rel": float(np.corrcoef(rel.ploidy, rel.median_abs_over_rel)[0, 1]),
        "ratio_median_abs_over_rel_to_ploidy": q(rel.median_abs_over_rel / rel.ploidy),
        "corr_ploidy_vs_mean_rel": float(np.corrcoef(rel.ploidy, rel.mean_rel)[0, 1]),
    }

    # sanity: known amplified / deleted genes
    name_of = dict(zip(m.ModelID, m.model_StrippedCellLineName))
    ln = {v: k for k, v in name_of.items()}
    checks = {"ERBB2 (2064)": ["SKBR3", "BT474", "HCC1954", "NCIN87", "OE19", "MCF7"],
              "MYC (4609)": ["COLO320", "HL60", "NCIH2171", "SW480", "HCT116"],
              "EGFR (1956)": ["A431", "MDAMB468", "HCC827", "A549"],
              "MDM2 (4193)": ["SJSA1", "NCIH460", "MCF7"],
              "CCNE1 (898)": ["OVCAR3", "HCC1569", "MCF7"],
              "CDKN2A (1029)": ["A549", "U87MG", "MCF7", "HCT116"]}
    san = {}
    for g, lines in checks.items():
        j = cn_cols.index(g)
        col = X[:, j]
        pct = stats.rankdata(np.nan_to_num(col, nan=np.nanmedian(col))) / len(col)
        rec = {"cohort_median": float(np.nanmedian(col)), "cohort_frac_gt_2": float(np.nanmean(col > 2)),
               "cohort_frac_lt_0.25": float(np.nanmean(col < 0.25))}
        ej = ex_cols.index(g) if g in ex_cols else None
        if ej is not None:
            ok = ~np.isnan(col) & ~np.isnan(E[:, ej])
            rec["cis_spearman"] = float(stats.spearmanr(col[ok], E[ok, ej]).statistic)
        for line in lines:
            if line in ln:
                i = m.ModelID.tolist().index(ln[line])
                rec[line] = {"rel_cn": None if np.isnan(col[i]) else round(float(col[i]), 3),
                             "percentile": round(float(pct[i]), 4)}
            else:
                rec[line] = "not matched"
        san[g] = rec
    R["sanity_known_cn"] = san

    # cis CN-expression correlation (all shared genes)
    bi_cn = [cn_cols.index(g) for g in both]
    bi_ex = [ex_cols.index(g) for g in both]
    Xc, Ec = X[:, bi_cn], E[:, bi_ex]
    rho = np.full(len(both), np.nan)
    for k in range(len(both)):
        ok = ~np.isnan(Xc[:, k]) & ~np.isnan(Ec[:, k])
        if ok.sum() > 50 and np.nanstd(Ec[ok, k]) > 0 and np.nanstd(Xc[ok, k]) > 0:
            rho[k] = stats.spearmanr(Xc[ok, k], Ec[ok, k]).statistic
    R["cis_spearman"] = {"quantiles": q(rho), "frac_positive": float(np.nanmean(rho > 0)),
                         "frac_gt_0.3": float(np.nanmean(rho > 0.3)), "n_genes": int(np.sum(~np.isnan(rho)))}
    pd.DataFrame({"depmap_column": both, "cis_spearman": rho}).to_csv(PROC / "depmap_cis_spearman.tsv", sep="\t", index=False)

    # ---------------- missingness ----------------
    gnan = np.isnan(X).mean(0)
    mnan = np.isnan(X).mean(1)
    R["missingness"] = {
        "cn_frac_missing_total": float(np.isnan(X).mean()),
        "cn_genes_all_missing": int((gnan == 1).sum()), "cn_genes_any_missing": int((gnan > 0).sum()),
        "cn_genes_gt_10pct_missing": int((gnan > 0.1).sum()), "cn_genes_gt_50pct_missing": int((gnan > 0.5).sum()),
        "cn_per_gene_missing_q": q(gnan), "cn_per_model_missing_q": q(mnan),
        "cn_models_gt_5pct_missing": int((mnan > 0.05).sum()),
        "cn_missing_by_chrom_class": {k: float(np.isnan(X[:, [cn_cols.index(c) for c in a.loc[cn_cols][a.loc[cn_cols].chrom_class == k].index]]).mean())
                                      for k in a.loc[cn_cols].chrom_class.dropna().unique()},
        "expr_frac_missing_total": float(np.isnan(E).mean()),
        "expr_genes_any_missing": int(np.isnan(E).any(0).sum()),
        "expr_genes_zero_in_all": int((np.nan_to_num(E) == 0).all(0).sum()),
        "expr_genes_median_lt_1": int((np.nanmedian(E, 0) < 1).sum()),
    }
    # Y-chromosome CN by sex (sanity)
    ycols = [c for c in cn_cols if a.loc[c, "chrom_class"] == "Y"]
    if ycols:
        yi = [cn_cols.index(c) for c in ycols]
        ymed = np.nanmedian(X[:, yi], 1)
        R["missingness"]["y_cn_median_by_sex"] = {s: q(ymed[(m.model_Sex == s).to_numpy()], (0.25, 0.5, 0.75))
                                                  for s in ["Male", "Female"]}
        R["missingness"]["y_cn_frac_missing_by_sex"] = {s: float(np.isnan(X[(m.model_Sex == s).to_numpy()][:, yi]).mean())
                                                        for s in ["Male", "Female"]}

    # ---------------- CN variability ----------------
    sd = np.nanstd(L, 0)
    mad = np.nanmedian(np.abs(L - np.nanmedian(L, 0)), 0)
    n03 = np.sum(np.abs(L) > 0.3, 0)
    n1 = np.sum(np.abs(L) > 1, 0)
    ngain = np.sum(L > 0.3, 0)
    nloss = np.sum(L < -0.3, 0)
    var_df = pd.DataFrame({"depmap_column": cn_cols, "sd_log2": sd, "mad_log2": mad, "n_abs_gt_0.3": n03,
                           "n_abs_gt_1": n1, "n_gain_0.3": ngain, "n_loss_0.3": nloss,
                           "frac_missing": gnan})
    var_df.to_csv(PROC / "depmap_cn_variability.tsv", sep="\t", index=False)
    N = len(m)
    defs = {
        "sd_log2>0.1": sd > 0.1, "sd_log2>0.2": sd > 0.2, "sd_log2>0.3": sd > 0.3, "sd_log2>0.5": sd > 0.5,
        "n|log2|>0.3>=10": n03 >= 10, "n|log2|>0.3>=20": n03 >= 20, "n|log2|>0.3>=50": n03 >= 50,
        "n|log2|>0.3>=5%": n03 >= 0.05 * N, "n|log2|>0.3>=10%": n03 >= 0.10 * N,
        "n|log2|>1>=10": n1 >= 10, "n|log2|>1>=20": n1 >= 20,
        "gain>=5%&loss>=5%": (ngain >= 0.05 * N) & (nloss >= 0.05 * N),
    }
    tg = {}
    for cell in ["hct116", "hek293t"]:
        cols = pie[pie[f"target_{cell}"] & pie.depmap_cn_column.notna()].depmap_cn_column.str.split(";").explode()
        tg[cell] = set(cols)
    R["cn_variability"] = {
        "sd_log2_q": q(sd), "mad_log2_q": q(mad), "n_abs_gt_0.3_q": q(n03), "n_abs_gt_1_q": q(n1),
        "n_models": N,
        "definitions": {k: {"all_cn_genes": int(msk.sum()),
                            "autosomal": int((msk & (a.loc[cn_cols].chrom_class == "autosome").to_numpy()).sum()),
                            **{f"pie_targets_{c}": int(sum(1 for g, mm in zip(cn_cols, msk) if mm and g in tg[c]))
                               for c in tg}}
                        for k, msk in defs.items()},
        "pie_targets_with_cn": {c: len(tg[c]) for c in tg},
        "note": "log2 of relative CN (floor 2^-5); relative CN is ploidy-normalised so |log2|>0.3 = gain/loss vs the line's own ploidy",
    }

    # ---------------- lineage / media / covariates ----------------
    lin = m.model_OncotreeLineage.value_counts()
    R["lineage"] = {"n_lineages": int(lin.size), "counts": {k: int(v) for k, v in lin.items()},
                    "n_lineages_ge_10": int((lin >= 10).sum()), "n_lineages_ge_20": int((lin >= 20).sum()),
                    "n_primary_disease": int(m.model_OncotreePrimaryDisease.nunique()),
                    "n_subtype": int(m.model_OncotreeSubtype.nunique()),
                    "n_oncotree_code": int(m.model_OncotreeCode.nunique()),
                    "missing_lineage": int(m.model_OncotreeLineage.isna().sum())}
    media = {}
    for c in ["model_OnboardedMedia", "model_FormulationID", "model_SerumFreeMedia", "model_GrowthPattern",
              "model_PlateCoating", "mc_DataSource", "mc_CellFormat", "mc_GrowthMedia", "mc_FormulationID",
              "mc_SerumFreeMedia", "mc_PlateCoating", "rna_stranded"]:
        s = m[c]
        cnt = s.value_counts()
        media[c] = {"n_levels": int(s.nunique()), "n_missing": int(s.isna().sum()),
                    "n_levels_ge_10": int((cnt >= 10).sum()),
                    "frac_in_levels_ge_10": float(cnt[cnt >= 10].sum() / len(s)),
                    "top": {str(k): int(v) for k, v in cnt.head(6).items()}}
    R["media"] = media
    cov = {}
    for c in ["sig24_MSIScore", "sig24_Ploidy", "sig24_WGD", "sig24_CIN", "sig24_LoHFraction", "sig24_Aneuploidy",
              "cmp_msi_status", "cmp_ploidy", "cmp_mutational_burden"]:
        cov[c] = {"n": int(m[c].notna().sum()), "frac": float(m[c].notna().mean()),
                  "by_model_type": {k: int(v) for k, v in m[m[c].notna()].model_ModelType.value_counts().items()}}
    cov["sig24_source"] = vc(m.sig24_source)
    cov["sig24_ploidy_q"] = q(m.sig24_Ploidy)
    cov["sig24_WGD_rate"] = float(m.sig24_WGD.mean())
    cov["sig24_MSI_gt20"] = int((m.sig24_MSIScore > 20).sum())
    ct = pd.crosstab(m.sig24_MSIScore > 20, m.cmp_msi_status).to_dict() if m.cmp_msi_status.notna().any() else {}
    sub = m[m.sig24_MSIScore.notna() & m.cmp_msi_status.notna()]
    cov["msi_concordance"] = {"n": len(sub), "crosstab(sig24>20 rows; cmp cols)":
                              {str(k): {str(kk): int(vv) for kk, vv in v.items()} for k, v in
                               pd.crosstab(sub.sig24_MSIScore > 20, sub.cmp_msi_status).to_dict().items()}}
    sub = m[m.sig24_Ploidy.notna() & m.cmp_ploidy.notna()]
    cov["ploidy_concordance"] = {"n": len(sub), "pearson": float(np.corrcoef(sub.sig24_Ploidy, sub.cmp_ploidy)[0, 1])}
    cov["any_ploidy"] = int((m.sig24_Ploidy.notna() | m.cmp_ploidy.notna()).sum())
    cov["any_msi"] = int((m.sig24_MSIScore.notna() | m.cmp_msi_status.notna()).sum())
    R["covariates"] = cov

    # ---------------- derived arm-level CN ----------------
    arm = pd.read_parquet(PROC / "depmap_arm_cn_derived.parquet")
    arm = arm[arm.ModelID.isin(m.ModelID)]
    W = arm.pivot(index="ModelID", columns="arm_id", values="median_log2_cn")
    ng = arm.groupby("arm_id").n_genes.median()
    W = W[[c for c in W.columns if ng[c] >= 20]]
    C = W.corr()
    chrom_of = {c: c[:-1] for c in C.columns}
    pairs = []
    for i, x in enumerate(C.columns):
        for y in C.columns[i + 1:]:
            if chrom_of[x] != chrom_of[y]:
                pairs.append((x, y, float(C.loc[x, y])))
    pairs.sort(key=lambda t: -abs(t[2]))
    R["arm_cn"] = {"n_arms_used": int(W.shape[1]), "arms_dropped_lt20_genes": [c for c in ng.index if ng[c] < 20],
                   "median_genes_per_arm": float(ng.median()),
                   "frac_model_arm_gain_gt_0.3": float((W > 0.3).mean().mean()),
                   "frac_model_arm_loss_lt_-0.3": float((W < -0.3).mean().mean()),
                   "interchrom_abs_corr_q": q([abs(p[2]) for p in pairs]),
                   "top_interchrom_pairs": [list(p) for p in pairs[:12]],
                   "n_interchrom_pairs_abs_r_gt_0.3": int(sum(abs(p[2]) > 0.3 for p in pairs))}
    # ploidy confounding: correlation of derived arms with 24Q4 ploidy (relative CN should be ~uncorrelated)
    pl_ = m.set_index("ModelID").sig24_Ploidy.reindex(W.index)
    ok = pl_.notna()
    R["arm_cn"]["max_abs_corr_arm_vs_ploidy"] = float(W[ok].corrwith(pl_[ok]).abs().max())

    out = PROC / "depmap_profile.json"
    with open(out, "w") as fh:
        json.dump(R, fh, indent=1, default=str)
    print(json.dumps(R, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
