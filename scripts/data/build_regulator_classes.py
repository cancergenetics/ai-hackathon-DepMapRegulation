#!/usr/bin/env python
"""Regulator classes from HGNC gene groups (D8 report analysis; rules in docs/BENCHMARK.md §11.1, fixed
before any result). Label-free: uses only HGNC and gene_annotation.tsv.

Output: data/processed/regulator_classes.tsv (gene key, hgnc_id, reg_class, matched groups) for every
DepMap gene key.
Usage: .venv/bin/python scripts/data/build_regulator_classes.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ORDER = ("ribosomal", "RBP", "TF", "kinase")
RULES = {
    "ribosomal": [re.compile(r"ribosomal proteins$")],
    "RBP": [re.compile(r"RNA binding"), re.compile(r"^Heterogeneous nuclear ribonucleoproteins$"),
            re.compile(r"^DEAD-box helicases$"), re.compile(r"^DEAH-box helicases$")],
    "TF": [re.compile(r"^Zinc fingers C2H2-type$"), re.compile(r"homeobox", re.I),
           re.compile(r"^Basic helix-loop-helix proteins$"), re.compile(r"bHLH transcription factor"),
           re.compile(r"^Basic leucine zipper proteins$"), re.compile(r"bZIP transcription factor"),
           re.compile(r"^Forkhead boxes$"), re.compile(r"^Nuclear receptor subfamily"),
           re.compile(r"^GATA zinc finger domain containing$")],
}
TF_GENERIC = re.compile(r"transcription factor", re.I)
TF_EXCLUDE = re.compile(r"^General transcription factor")
KINASE = re.compile(r"kinase", re.I)
KINASE_EXCLUDE = re.compile(r"anchoring|kinase phosphatase|kinase inhibitor|regulatory subunit|"
                            r"kinase activators|non-catalytic", re.I)


def group_classes(group: str) -> set:
    out = {c for c, rs in RULES.items() if any(r.search(group) for r in rs)}
    if TF_GENERIC.search(group) and not TF_EXCLUDE.search(group):
        out.add("TF")
    if KINASE.search(group) and not KINASE_EXCLUDE.search(group):
        out.add("kinase")
    return out


def classify(groups: str | float) -> tuple[str, str]:
    if not isinstance(groups, str) or not groups.strip():
        return "other", ""
    gs = [g.strip() for g in groups.split("|")]
    hit = {g: group_classes(g) for g in gs}
    cls = set().union(*hit.values())
    for c in ORDER:
        if c in cls:
            return c, "|".join(g for g, h in hit.items() if c in h)
    return "other", ""


def main() -> int:
    hg = pd.read_csv(ROOT / "starting_data" / "gene_with_protein_product.txt", sep="\t", low_memory=False,
                     usecols=["hgnc_id", "symbol", "gene_group"])
    ann = pd.read_csv(ROOT / "data" / "processed" / "gene_annotation.tsv", sep="\t", low_memory=False,
                      usecols=["depmap_column", "hgnc_id"])
    m = ann.merge(hg[["hgnc_id", "gene_group"]], on="hgnc_id", how="left")
    if m["depmap_column"].duplicated().any():
        raise SystemExit("HGNC id maps to several HGNC rows")
    res = m["gene_group"].map(classify)
    out = pd.DataFrame({"gene": m["depmap_column"], "hgnc_id": m["hgnc_id"],
                        "reg_class": [r[0] for r in res], "matched_groups": [r[1] for r in res]})
    p = ROOT / "data" / "processed" / "regulator_classes.tsv"
    out.to_csv(p, sep="\t", index=False)
    print(out["reg_class"].value_counts().to_string(), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
