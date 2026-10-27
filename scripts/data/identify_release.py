#!/usr/bin/env python
"""Collect evidence for which DepMap release starting_data/ comes from. Writes
data/processed/release_evidence.json.

The authoritative check would compare md5s with the DepMap portal file listing
(https://depmap.org/portal/api/download/files), but that endpoint sits behind a human-verification
page and cannot be scripted. Evidence therefore comes from (a) figshare deposits, which are md5-comparable
but end at 24Q4, and (b) release-note markers checked in the files themselves:
  - 24Q4 figshare md5 of Model.csv / ModelCondition.csv / OmicsProfiles.csv  (mismatch => not 24Q4)
  - metadata columns (SequencingID, ModelConditionID, IsDefaultEntryForMC) in omics matrices (25Q3+)
  - SubtypeTree.csv present (25Q3+); PAM50 breast subtype nodes ("new breast cancer subtypes", 26Q3)
  - latest SequencingDate in OmicsProfiles.csv
  - presence of omics for the NextGen 3D models (Neiswender, Maffa, Brenan et al., Nature 2026),
    integrated into the portal in 26Q3 (figshare 29472362 model_metadata.csv, CC BY 4.0)
"""
from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import pandas as pd
import truststore

truststore.inject_into_ssl()
ROOT = Path(__file__).resolve().parents[2]
SD = ROOT / "starting_data"
NG_URL = "https://ndownloader.figshare.com/files/64494438"  # figshare 29472362 model_metadata.csv
NG_MD5 = "97b26ec52b8a05e52df7abd507b8139f"
NG = ROOT / "data" / "external" / "nextgen_2026" / "model_metadata.csv"
FIGSHARE_24Q4 = 27993248


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        while c := fh.read(1 << 22):
            h.update(c)
    return h.hexdigest()


def main() -> int:
    ev: dict = {}
    with urllib.request.urlopen(f"https://api.figshare.com/v2/articles/{FIGSHARE_24Q4}/files?page_size=500") as r:
        f24 = {f["name"]: f["supplied_md5"] for f in json.load(r)}
    ev["md5_equal_to_24Q4_figshare"] = {n: md5(SD / n) == f24.get(n) for n in ["Model.csv", "ModelCondition.csv", "OmicsProfiles.csv"]}
    hdr = pd.read_csv(SD / "OmicsCNGeneWGS.csv", nrows=0).columns.tolist()
    ev["omics_matrix_metadata_columns"] = [c for c in hdr[:7] if not c.endswith(")")]
    st = pd.read_csv(SD / "SubtypeTree.csv")
    ev["subtype_tree_pam50_nodes"] = sorted(st.NodeName[st.NodeName.str.contains("PAM50", na=False)].tolist())
    prof = pd.read_csv(SD / "OmicsProfiles.csv")
    d = pd.to_datetime(prof.SequencingDate, errors="coerce")
    ev["max_sequencing_date_by_datatype"] = {k: str(v.date()) for k, v in d.groupby(prof.DataType).max().items()}
    if not NG.exists():
        NG.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(NG_URL, NG)
    assert md5(NG) == NG_MD5
    ng = pd.read_csv(NG)
    ng = ng[ng.IsNextGen == True]  # noqa: E712
    wgs = set(prof[prof.DataType == "wgs"].ModelID)
    rna = set(prof[prof.DataType == "rna"].ModelID)
    ev["nextgen_models"] = {"n": len(ng), "with_HasWGS": int(ng.HasWGS.sum()),
                            "HasWGS_with_wgs_profile_here": int(ng[ng.HasWGS == True].ModelID.isin(wgs).sum()),  # noqa: E712
                            "with_rna_profile_here": int(ng.ModelID.isin(rna).sum())}
    model = pd.read_csv(SD / "Model.csv")
    ev["model_csv"] = {"n_models": len(model), "model_type": model.ModelType.value_counts().to_dict()}
    ev["conclusion"] = ("DepMap Public 26Q3 (portal release announced 2026-10-05): not 24Q4 (md5 differ); has 25Q3+ "
                        "omics metadata columns and SubtypeTree; PAM50 breast subtypes (26Q3 release note); RNA "
                        "sequenced up to Oct 2025; NextGen-manuscript models present. Not md5-verified against the "
                        "portal listing (bot-protected).")
    out = ROOT / "data" / "processed" / "release_evidence.json"
    out.write_text(json.dumps(ev, indent=1, default=str))
    print(json.dumps(ev, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
