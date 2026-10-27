#!/usr/bin/env python
"""Human paralogue pairs from Ensembl release 116 (the release of gene_annotation.tsv), for the
paralog-pair exclusion sensitivity analysis (DESIGN §5.5 / §13; D8; docs/BENCHMARK.md §11.2).

Source: the Ensembl Compara 116 homology TSV dump for Homo sapiens,
  https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/homo_sapiens/
  Compara.116.protein_default.homologies.tsv.gz
It is verified against the Ensembl MD5SUM file. This is the gene-tree homology set that BioMart's
`hsapiens_paralog_ensembl_gene` attribute serves. The BioMart query (jun2026.archive.ensembl.org, Ensembl
Genes 116) did not complete in about 40 min on 2026-10-27, so the release-pinned, checksummed dump is used.
Paralogues: rows with homology_species == homo_sapiens and homology_type in {within_species_paralog,
other_paralog} (gene_split excluded).
Licence: Ensembl, no restrictions on use.

Outputs:
  data/external/annotation/Compara.116.protein_default.homologies.homo_sapiens.tsv.gz  (raw download)
  data/processed/paralog_pairs.tsv  symmetric gene-key pairs (gene_a, gene_b, homology_type, identity)
  a row in data/external/annotation_manifest.tsv
Usage: .venv/bin/python scripts/data/download_paralogs.py [--skip-download]
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import sys
import urllib.request
from pathlib import Path

import pandas as pd
import truststore

truststore.inject_into_ssl()      # system trust store (as in download_annotation.py)

ROOT = Path(__file__).resolve().parents[2]
RELEASE = "116"
BASE = f"https://ftp.ensembl.org/pub/release-{RELEASE}/tsv/ensembl-compara/homologies/homo_sapiens/"
FILE = f"Compara.{RELEASE}.protein_default.homologies.tsv.gz"
RAW = ROOT / "data" / "external" / "annotation" / f"Compara.{RELEASE}.protein_default.homologies.homo_sapiens.tsv.gz"
OUT = ROOT / "data" / "processed" / "paralog_pairs.tsv"
MANIFEST = ROOT / "data" / "external" / "annotation_manifest.tsv"
PARALOG_TYPES = {"within_species_paralog", "other_paralog"}


def digest(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def download() -> None:
    md5sum = urllib.request.urlopen(BASE + "MD5SUM", timeout=120).read().decode()
    want = {ln.split()[1]: ln.split()[0] for ln in md5sum.splitlines() if ln.strip()}[FILE]
    tmp = RAW.with_suffix(".part")
    with urllib.request.urlopen(BASE + FILE, timeout=1800) as r, open(tmp, "wb") as out:
        while True:
            b = r.read(1 << 22)
            if not b:
                break
            out.write(b)
    if digest(tmp, "md5") != want:
        tmp.unlink()
        raise SystemExit(f"md5 mismatch for {FILE}")
    tmp.rename(RAW)
    row = ["ensembl_compara_homologies", RAW.name, BASE + FILE, f"Ensembl Compara release {RELEASE} (protein_default)",
           str(RAW.stat().st_size), digest(RAW, "sha256"), want, "md5 matches Ensembl MD5SUM",
           dt.date.today().isoformat()]
    with open(MANIFEST, "a") as fh:
        fh.write("\t".join(row) + "\n")


def build() -> None:
    cols = ["gene_stable_id", "species", "identity", "homology_type", "homology_gene_stable_id", "homology_species"]
    parts = []
    for ch in pd.read_csv(RAW, sep="\t", usecols=cols, chunksize=2_000_000):
        ch = ch[(ch["species"] == "homo_sapiens") & (ch["homology_species"] == "homo_sapiens")
                & ch["homology_type"].isin(PARALOG_TYPES)]
        parts.append(ch)
    h = pd.concat(parts, ignore_index=True)
    ann = pd.read_csv(ROOT / "data" / "processed" / "gene_annotation.tsv", sep="\t", low_memory=False,
                      usecols=["depmap_column", "ensembl_gene_id"]).dropna()
    if ann["ensembl_gene_id"].duplicated().any():
        ann = ann.drop_duplicates("ensembl_gene_id", keep=False)    # ambiguous Ensembl ids are not used
    key = dict(zip(ann["ensembl_gene_id"], ann["depmap_column"]))
    a = h["gene_stable_id"].map(key)
    b = h["homology_gene_stable_id"].map(key)
    ok = a.notna() & b.notna() & (a != b)
    d = pd.DataFrame({"gene_a": a[ok], "gene_b": b[ok], "homology_type": h["homology_type"][ok],
                      "identity": h["identity"][ok]})
    d = pd.concat([d, d.rename(columns={"gene_a": "gene_b", "gene_b": "gene_a"})])
    d = d.sort_values(["gene_a", "gene_b", "identity"], ascending=[True, True, False]).drop_duplicates(["gene_a", "gene_b"])
    d.to_csv(OUT, sep="\t", index=False)
    print(f"{len(h)} human paralogue rows -> {len(d)} directed DepMap gene-key pairs "
          f"({d['gene_a'].nunique()} genes); {OUT}", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-download", action="store_true")
    a = ap.parse_args()
    if not a.skip_download:
        download()
    build()
    return 0


if __name__ == "__main__":
    sys.exit(main())
