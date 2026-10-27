#!/usr/bin/env python
"""Download pinned gene-annotation resources and record URLs + checksums.

- HGNC complete set, monthly archive (pinned date) -- CC0. Needed in addition to
  starting_data/gene_with_protein_product.txt because Perturb-seq features include non-coding genes.
- Ensembl GRCh38 GTF, pinned release (gene coordinates, biotypes, stable IDs). Verified against the
  Ensembl CHECKSUMS file (BSD `sum` algorithm).
- UCSC hg38 cytoBand.txt.gz (cytoband + centromere (acen) positions for p/q arm assignment).

Usage: .venv/bin/python scripts/data/download_annotation.py
Writes data/external/annotation/* and data/external/annotation_manifest.tsv
"""
from __future__ import annotations

import base64
import datetime as dt
import hashlib
import sys
import urllib.request
from pathlib import Path

import truststore

truststore.inject_into_ssl()

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "external" / "annotation"
MANIFEST = ROOT / "data" / "external" / "annotation_manifest.tsv"

HGNC_DATE = "2026-10-06"
ENSEMBL_RELEASE = 116
HGNC_URL = ("https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/"
            f"hgnc_complete_set_{HGNC_DATE}.tsv")
HGNC_META = ("https://storage.googleapis.com/storage/v1/b/public-download-files/o/"
             f"hgnc%2Farchive%2Farchive%2Fmonthly%2Ftsv%2Fhgnc_complete_set_{HGNC_DATE}.tsv?fields=md5Hash,size")
ENS_BASE = f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/gtf/homo_sapiens/"
GTF_NAME = f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz"
CYTO_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as out:
        while chunk := r.read(1 << 22):
            out.write(chunk)
    tmp.rename(dest)


def read_url(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def hashes(path: Path) -> tuple[str, str, int]:
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 22):
            sha.update(chunk)
            md5.update(chunk)
    return sha.hexdigest(), md5.hexdigest(), path.stat().st_size


def bsd_sum(path: Path) -> tuple[int, int]:
    """BSD `sum` checksum (as used in Ensembl CHECKSUMS files): (checksum, 1 KiB blocks)."""
    s = 0
    n = 0
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            n += len(chunk)
            for b in chunk:
                s = (s >> 1) + ((s & 1) << 15)
                s = (s + b) & 0xFFFF
    return s, (n + 1023) // 1024


def main() -> int:
    import json

    rows = []
    today = dt.date.today().isoformat()

    # HGNC
    dest = OUT / f"hgnc_complete_set_{HGNC_DATE}.tsv"
    if not dest.exists():
        fetch(HGNC_URL, dest)
    meta = json.loads(read_url(HGNC_META))
    sha, md5, size = hashes(dest)
    gcs_md5 = base64.b64decode(meta["md5Hash"]).hex()
    if gcs_md5 != md5:
        raise SystemExit(f"HGNC md5 mismatch {md5} vs {gcs_md5}")
    rows.append(("hgnc", dest.name, HGNC_URL, f"HGNC monthly {HGNC_DATE}", size, sha, md5,
                 "md5 matches GCS object md5", today))

    # Ensembl GTF
    dest = OUT / GTF_NAME
    if not dest.exists():
        print(f"downloading {GTF_NAME}", file=sys.stderr)
        fetch(ENS_BASE + GTF_NAME, dest)
    checks = read_url(ENS_BASE + "CHECKSUMS").decode().splitlines()
    want = next(l.split() for l in checks if l.split()[-1] == GTF_NAME)
    got = bsd_sum(dest)
    if (int(want[0]), int(want[1])) != got:
        raise SystemExit(f"Ensembl GTF BSD sum mismatch {got} vs {want}")
    sha, md5, size = hashes(dest)
    rows.append(("ensembl", GTF_NAME, ENS_BASE + GTF_NAME, f"Ensembl release {ENSEMBL_RELEASE} (GRCh38)",
                 size, sha, md5, "BSD sum matches Ensembl CHECKSUMS", today))

    # UCSC cytoBand
    dest = OUT / "cytoBand.hg38.txt.gz"
    if not dest.exists():
        fetch(CYTO_URL, dest)
    sha, md5, size = hashes(dest)
    rows.append(("ucsc", dest.name, CYTO_URL, "hg38 cytoBand (Last-Modified recorded by server)", size, sha,
                 md5, "no upstream checksum", today))

    with open(MANIFEST, "w") as fh:
        fh.write("source\tfile\turl\tversion\tsize_bytes\tsha256\tmd5\tverification\tretrieved\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {MANIFEST}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
