#!/usr/bin/env python
"""Build analysis-ready DepMap tables (no modelling).

Matching rule: for each modality take the profile row flagged IsDefaultEntryForModel == "Yes"
(exactly one per model in 26Q3 for both CN and RNA; asserted). A model is "matched" if it has a
default CN (WGS) row and a default RNA row. CN and RNA default rows may come from different
ModelConditions (flagged `same_mc`).

Outputs (data/processed/):
  depmap_models.tsv           one row per model with a default CN and/or RNA profile: profile ids,
                              Model.csv / ModelCondition.csv annotations, 24Q4 PureCN/MSIsensor2
                              signatures (with provenance), Sanger Cell Model Passports MSI/ploidy.
  depmap_cn_matched.parquet   ModelID + all CN gene columns, matched models only, sorted by ModelID
  depmap_expr_matched.parquet ModelID + all expression gene columns, same row order
  depmap_arm_cn_derived.parquet  long table (ModelID, chrom, arm, arm_id, n_genes, median_log2_cn,
                              mean_log2_cn) for all default-CN models; derived from gene-level
                              relative CN because 26Q3 ships no arm-level calls for WGS profiles.
                              log2 uses a floor of CN_FLOOR (relative CN) to keep deletions finite.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
SD = ROOT / "starting_data"
EXT = ROOT / "data" / "external"
PROC = ROOT / "data" / "processed"
INTERIM = ROOT / "data" / "interim"
META = ["SequencingID", "ModelConditionID", "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel"]
CN_FLOOR = 2 ** -5
MODEL_COLS = ["CellLineName", "StrippedCellLineName", "DepmapModelType", "ModelType", "OncotreeLineage",
              "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode", "Sex", "Age", "AgeCategory",
              "PrimaryOrMetastasis", "GrowthPattern", "OnboardedMedia", "FormulationID", "SerumFreeMedia",
              "PlateCoating", "EngineeredModel", "PediatricModelType", "RRID", "SangerModelID", "CCLEName"]
MC_COLS = ["DataSource", "CellFormat", "GrowthMedia", "FormulationID", "SerumFreeMedia", "PlateCoating"]
SIGS = ["MSIScore", "Ploidy", "LoHFraction", "WGD", "CIN", "Aneuploidy"]


def default_rows(path: Path) -> pd.DataFrame:
    m = pl.read_parquet(path, columns=META).to_pandas()
    d = m[m.IsDefaultEntryForModel == "Yes"]
    assert not d.ModelID.duplicated().any(), f"multiple default rows per model in {path.name}"
    return d


def main() -> int:
    prof = pd.read_csv(SD / "OmicsProfiles.csv")
    seq2prof = prof.set_index("SequencingID")
    cn_d = default_rows(INTERIM / "depmap_cn_wgs.parquet")
    rna_d = default_rows(INTERIM / "depmap_expr_tpm_logp1.parquet")
    cn_d = cn_d.rename(columns={"SequencingID": "cn_seq_id", "ModelConditionID": "cn_mc"})[["ModelID", "cn_seq_id", "cn_mc"]]
    rna_d = rna_d.rename(columns={"SequencingID": "rna_seq_id", "ModelConditionID": "rna_mc"})[["ModelID", "rna_seq_id", "rna_mc"]]
    t = cn_d.merge(rna_d, on="ModelID", how="outer")
    t["has_cn"] = t.cn_seq_id.notna()
    t["has_rna"] = t.rna_seq_id.notna()
    t["matched"] = t.has_cn & t.has_rna
    t["same_mc"] = np.where(t.matched, t.cn_mc == t.rna_mc, pd.NA)
    t["cn_profile_id"] = t.cn_seq_id.map(seq2prof.ProfileID)
    t["rna_profile_id"] = t.rna_seq_id.map(seq2prof.ProfileID)
    t["cn_seq_date"] = t.cn_seq_id.map(seq2prof.SequencingDate)
    t["rna_seq_date"] = t.rna_seq_id.map(seq2prof.SequencingDate)
    t["rna_stranded"] = t.rna_seq_id.map(seq2prof.Stranded)

    model = pd.read_csv(SD / "Model.csv").set_index("ModelID")
    for c in MODEL_COLS:
        t[f"model_{c}"] = t.ModelID.map(model[c])
    mc = pd.read_csv(SD / "ModelCondition.csv").set_index("ModelConditionID")
    # media/culture annotation of the condition the RNA was profiled in (expression is the response)
    mc_key = t.rna_mc.fillna(t.cn_mc)
    for c in MC_COLS:
        t[f"mc_{c}"] = mc_key.map(mc[c])

    # 24Q4 signatures: prefer the value for the *same* profile as the 26Q3 default CN profile
    sp = pd.read_csv(EXT / "depmap_24Q4" / "OmicsSignaturesProfile.csv", index_col=0)
    sm = pd.read_csv(EXT / "depmap_24Q4" / "OmicsSignatures.csv", index_col=0)
    same = t.cn_profile_id.isin(sp.index)
    t["sig24_source"] = np.where(same, "same_profile", np.where(t.ModelID.isin(sm.index), "model_level_24Q4", "none"))
    for c in SIGS:
        v_prof = t.cn_profile_id.map(sp[c])
        v_mod = t.ModelID.map(sm[c])
        t[f"sig24_{c}"] = np.where(same, v_prof, v_mod)
    # 24Q4 profile-level datatype of the profile the model-level signature came from (WES vs WGS)
    p24 = pd.read_csv(EXT / "depmap_24Q4" / "OmicsProfiles.csv").set_index("ProfileID")
    d24 = pd.read_csv(EXT / "depmap_24Q4" / "OmicsDefaultModelProfiles.csv")
    d24 = d24[d24.ProfileType == "dna"].drop_duplicates("ModelID").set_index("ModelID")
    t["sig24_model_profile"] = t.ModelID.map(d24.ProfileID)
    t["sig24_model_profile_datatype"] = t.sig24_model_profile.map(p24.Datatype)

    # Sanger Cell Model Passports (independent annotation)
    cmp_ = pd.read_csv(EXT / "sanger_cmp" / "model_list_latest.csv.gz", low_memory=False)
    cmp_ = cmp_.drop_duplicates("model_id").set_index("model_id")
    t["cmp_msi_status"] = t.model_SangerModelID.map(cmp_.msi_status)
    t["cmp_ploidy"] = t.model_SangerModelID.map(cmp_.ploidy)
    t["cmp_mutational_burden"] = t.model_SangerModelID.map(cmp_.mutational_burden)

    t = t.sort_values("ModelID").reset_index(drop=True)
    t.to_csv(PROC / "depmap_models.tsv", sep="\t", index=False)

    # matched matrices
    keep = t[t.matched]
    for src, dest, key in [("depmap_cn_wgs.parquet", "depmap_cn_matched.parquet", "cn_seq_id"),
                           ("depmap_expr_tpm_logp1.parquet", "depmap_expr_matched.parquet", "rna_seq_id")]:
        lf = pl.scan_parquet(INTERIM / src)
        df = (lf.filter(pl.col("SequencingID").is_in(keep[key].tolist()))
              .drop([c for c in META if c not in ("ModelID",)]).collect())
        df = df.sort("ModelID")
        assert df["ModelID"].to_list() == keep.ModelID.tolist()
        df.write_parquet(PROC / dest, compression="zstd")
        print(f"{dest}: {df.shape}", file=sys.stderr)

    # derived arm-level CN for all default-CN models
    ann = pd.read_csv(PROC / "gene_annotation.tsv", sep="\t", dtype={"chrom": str})
    ann = ann[ann.in_cn & ann.arm.isin(["p", "q"]) & ann.chrom.notna() & (ann.chrom != "MT") & ~ann.is_par]
    cn = pl.read_parquet(INTERIM / "depmap_cn_wgs.parquet")
    cn = cn.filter(pl.col("SequencingID").is_in(cn_d.cn_seq_id.tolist())).sort("ModelID")
    models = cn["ModelID"].to_list()
    X = cn.select(ann.depmap_column.tolist()).to_numpy().astype(np.float64)
    L = np.log2(np.maximum(X, CN_FLOOR))
    out = []
    ann = ann.reset_index(drop=True)
    for (chrom, arm), idx in ann.groupby(["chrom", "arm"]).groups.items():
        sub = L[:, np.asarray(idx)]
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                med = np.nanmedian(sub, axis=1)
                mean = np.nanmean(sub, axis=1)
        n = np.sum(~np.isnan(sub), axis=1)
        out.append(pd.DataFrame({"ModelID": models, "chrom": chrom, "arm": arm, "arm_id": f"{chrom}{arm}",
                                 "n_genes": n, "median_log2_cn": med, "mean_log2_cn": mean}))
    arm = pd.concat(out, ignore_index=True)
    arm.to_parquet(PROC / "depmap_arm_cn_derived.parquet", index=False)
    print(f"depmap_models.tsv: {len(t)} models ({int(t.matched.sum())} matched); arms: {arm.arm_id.nunique()}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
