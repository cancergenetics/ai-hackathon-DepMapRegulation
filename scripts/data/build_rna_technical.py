#!/usr/bin/env python
"""RNA technical covariates for each model's default RNA profile (DESIGN v2 data task T4).

Join key: the default RNA SequencingID (unique in 26Q3 OmicsProfiles.csv). ProfileID is NOT unique there
(one ProfileID can carry several SequencingIDs = re-sequencing runs), so it is not used as a key.
Columns from 26Q3 OmicsProfiles.csv: Stranded, SequencingPlatform, SequencingDate, SourceModelCondition,
SharedToDbGaP, CellFormat, GrowthPattern. Column from 24Q4 OmicsProfiles.csv (matched by ProfileID, so
provenance-flagged): Product (sequencer) and Stranded (for a consistency check).
Output: data/processed/depmap_rna_technical.tsv (one row per model with a default RNA profile).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data" / "processed"


def main() -> int:
    models = pd.read_csv(PROC / "depmap_models.tsv", sep="\t", low_memory=False)
    prof = pd.read_csv(ROOT / "starting_data" / "OmicsProfiles.csv")
    assert prof.SequencingID.is_unique
    p = prof.set_index("SequencingID")
    t = models.loc[models.has_rna, ["ModelID", "rna_seq_id", "rna_profile_id", "rna_mc", "matched",
                                    "model_ModelType"]].copy()
    for c in ["Stranded", "SequencingPlatform", "SequencingDate", "SourceModelCondition", "SharedToDbGaP",
              "CellFormat", "GrowthPattern"]:
        t[f"rna_{c}"] = t.rna_seq_id.map(p[c])
    t["rna_seq_year"] = pd.to_datetime(t.rna_SequencingDate, errors="coerce").dt.year
    n_seq = prof.groupby("ProfileID").SequencingID.nunique()
    t["rna_profile_n_sequencings"] = t.rna_profile_id.map(n_seq)
    p24 = pd.read_csv(ROOT / "data" / "external" / "depmap_24Q4" / "OmicsProfiles.csv")
    p24 = p24[p24.Datatype == "rna"].drop_duplicates("ProfileID").set_index("ProfileID")
    t["rna24_Product"] = t.rna_profile_id.map(p24.Product)
    t["rna24_Stranded"] = t.rna_profile_id.map(p24.Stranded)
    t.to_csv(PROC / "depmap_rna_technical.tsv", sep="\t", index=False)
    both = t.rna24_Stranded.notna()
    print(f"{len(t)} models; Stranded missing {int(t.rna_Stranded.isna().sum())}; "
          f"24Q4 vs 26Q3 Stranded disagree {int((t.rna24_Stranded[both].astype(bool) != t.rna_Stranded[both].astype(bool)).sum())}/{int(both.sum())}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
