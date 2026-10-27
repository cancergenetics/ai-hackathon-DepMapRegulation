#!/usr/bin/env python
"""TP53 functional status per model from DepMap somatic mutation calls (DESIGN v2 data task T1).

Source: DepMap Public 24Q4 (figshare+ 27993248). This is the latest release whose mutation files can be
downloaded programmatically; 26Q3 mutation files are portal-only (bot-protected). Files:
  OmicsSomaticMutationsMatrixDamaging.csv  (LikelyLoF variants; 0 none, 1 AF sum <= 0.95, 2 AF sum > 0.95)
  OmicsSomaticMutationsMatrixHotspot.csv   (Hess 2019 / OncoKB / COSMIC tier-1 hotspots; same coding)
  OmicsSomaticMutations.csv                (MAF-like Mutect2 calls; used for detail columns + consistency)

Definitions (explicit):
  genotyped        := ModelID is a row of both 24Q4 matrices (models absent there have no 24Q4 calls)
  tp53_damaging    := damaging matrix value for "TP53 (7157)" (0/1/2), NA if not genotyped
  tp53_hotspot     := hotspot matrix value (0/1/2), NA if not genotyped
  tp53_status      := "mutant" if tp53_damaging > 0 or tp53_hotspot > 0; "wt" if genotyped and both 0;
                      "unknown" otherwise. (This is the DESIGN v2 covariate; a missingness indicator is
                      the modelling step's job.)
  tp53_other_coding := genotyped, status "wt", but the MAF has a TP53 coding variant that is neither
                      LikelyLoF nor hotspot (informational; not part of tp53_status)
  tp53_cn_homdel   := 26Q3 relative CN of TP53 < 0.25 (deep deletion; informational, from the 26Q3 CN)
Output: data/processed/depmap_tp53.tsv, one row per model in depmap_models.tsv.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
D24 = ROOT / "data" / "external" / "depmap_24Q4"
PROC = ROOT / "data" / "processed"
GENE = "TP53 (7157)"
CODING = ("missense", "stop_gained", "frameshift", "splice_acceptor", "splice_donor", "inframe", "start_lost",
          "stop_lost", "protein_altering")


def read_col(path: Path) -> pd.Series:
    first = pd.read_csv(path, nrows=0).columns[0]
    df = pd.read_csv(path, usecols=[first, GENE], index_col=0)
    assert df.index.is_unique
    return df[GENE]


def main() -> int:
    models = pd.read_csv(PROC / "depmap_models.tsv", sep="\t", low_memory=False)
    dam = read_col(D24 / "OmicsSomaticMutationsMatrixDamaging.csv")
    hot = read_col(D24 / "OmicsSomaticMutationsMatrixHotspot.csv")
    if set(dam.index) != set(hot.index):
        print(f"NOTE: damaging/hotspot matrices differ in rows ({len(dam)} vs {len(hot)})", file=sys.stderr)
    genotyped = set(dam.index) & set(hot.index)
    maf = pd.read_csv(D24 / "OmicsSomaticMutations.csv", low_memory=False,
                      usecols=["ModelID", "HugoSymbol", "ProteinChange", "VariantInfo", "LikelyLoF", "Hotspot", "AF"])
    maf = maf[maf.HugoSymbol == "TP53"]
    coding = maf.VariantInfo.fillna("").str.contains("|".join(CODING))
    det = maf.groupby("ModelID").agg(
        tp53_n_variants=("ProteinChange", "size"),
        tp53_protein_changes=("ProteinChange", lambda x: ";".join(sorted(set(map(str, x))))),
        tp53_max_af=("AF", "max"))
    other = maf[coding & ~maf.LikelyLoF.astype(bool) & ~maf.Hotspot.astype(bool)].ModelID.unique()

    t = models[["ModelID"]].copy()
    t["genotyped_24Q4"] = t.ModelID.isin(genotyped)
    t["tp53_damaging"] = t.ModelID.map(dam).where(t.genotyped_24Q4)
    t["tp53_hotspot"] = t.ModelID.map(hot).where(t.genotyped_24Q4)
    mut = (t.tp53_damaging > 0) | (t.tp53_hotspot > 0)
    t["tp53_status"] = np.where(~t.genotyped_24Q4, "unknown", np.where(mut, "mutant", "wt"))
    t = t.join(det, on="ModelID")
    t["tp53_n_variants"] = t.tp53_n_variants.fillna(0).astype(int).where(t.genotyped_24Q4)
    t["tp53_other_coding"] = t.genotyped_24Q4 & (t.tp53_status == "wt") & t.ModelID.isin(other)

    cn = pl.read_parquet(ROOT / "data" / "interim" / "depmap_cn_wgs.parquet",
                         columns=["SequencingID", "ModelID", "IsDefaultEntryForModel", GENE])
    cn = cn.filter(pl.col("IsDefaultEntryForModel") == "Yes").to_pandas().set_index("ModelID")[GENE]
    t["tp53_rel_cn_26Q3"] = t.ModelID.map(cn)
    t["tp53_cn_homdel"] = (t.tp53_rel_cn_26Q3 < 0.25).where(t.tp53_rel_cn_26Q3.notna())
    t["source"] = "DepMap Public 24Q4 OmicsSomaticMutationsMatrix{Damaging,Hotspot} (figshare+ 27993248)"

    # consistency: MAF LikelyLoF / Hotspot rows vs matrices
    lof_maf = set(maf[maf.LikelyLoF.astype(bool)].ModelID)
    hot_maf = set(maf[maf.Hotspot.astype(bool)].ModelID)
    g = t[t.genotyped_24Q4]
    print(f"genotyped {len(g)}; damaging>0 vs MAF LikelyLoF disagree: "
          f"{int(((g.tp53_damaging > 0) != g.ModelID.isin(lof_maf)).sum())}; hotspot disagree: "
          f"{int(((g.tp53_hotspot > 0) != g.ModelID.isin(hot_maf)).sum())}", file=sys.stderr)
    t.to_csv(PROC / "depmap_tp53.tsv", sep="\t", index=False)
    print(t.tp53_status.value_counts().to_dict(), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
