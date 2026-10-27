#!/usr/bin/env python
"""Record sha256 / md5 / size of every input dataset file in docs/data_checksums.tsv (committed).

Covers starting_data/*, the PIE X-Atlas/Orion files (sha256 also checked against the Hugging Face
LFS etag stored in data/external/pie_x_atlas_orion/.cache/huggingface/download/**/*.metadata), and
the files listed in data/external/{download,annotation}_manifest.tsv (written by the download scripts), and
data/external/nextgen_2026/* (identify_release.py).

Usage: .venv/bin/python scripts/data/checksum_inputs.py
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PIE = ROOT / "data" / "external" / "pie_x_atlas_orion"
OUT = ROOT / "docs" / "data_checksums.tsv"


def hashes(path: Path) -> tuple[str, str, int]:
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 22):
            sha.update(chunk)
            md5.update(chunk)
    return sha.hexdigest(), md5.hexdigest(), path.stat().st_size


def main() -> int:
    rows = []
    for p in sorted((ROOT / "starting_data").iterdir()):
        if p.is_file() and not p.name.startswith("."):
            sha, md5, size = hashes(p)
            rows.append(dict(dataset="depmap_26Q3_starting_data" if p.suffix == ".csv" else "hgnc_protein_coding",
                             path=str(p.relative_to(ROOT)), size_bytes=size, sha256=sha, md5=md5, verification=""))
    meta_dir = PIE / ".cache" / "huggingface" / "download"
    for p in sorted(PIE.rglob("*")):
        if not p.is_file() or ".cache" in p.parts:
            continue
        sha, md5, size = hashes(p)
        ver = ""
        mfile = meta_dir / p.relative_to(PIE).with_name(p.name + ".metadata")
        if mfile.exists():
            lines = mfile.read_text().split("\n")
            commit, etag = lines[0].strip(), lines[1].strip()
            if len(etag) == 64:
                ver = f"HF commit {commit}; sha256 {'==' if etag == sha else '!='} HF LFS etag"
            else:
                ver = f"HF commit {commit}; etag {etag} (git blob id, not sha256)"
        rows.append(dict(dataset="pie_x_atlas_orion", path=str(p.relative_to(ROOT)), size_bytes=size, sha256=sha,
                         md5=md5, verification=ver))
    for man in ["download_manifest.tsv", "annotation_manifest.tsv"]:
        mp = ROOT / "data" / "external" / man
        if not mp.exists():
            print(f"WARNING: {mp} missing; run the download scripts", file=sys.stderr)
            continue
        m = pd.read_csv(mp, sep="\t")
        for r in m.itertuples():
            sub = {"depmap_24Q4": "depmap_24Q4", "sanger_cmp": "sanger_cmp", "depmap_26Q3": "depmap_26Q3"}
            folder = sub.get(r.source, "annotation")
            p = ROOT / "data" / "external" / folder / r.file
            sha, md5, size = hashes(p)
            if sha != r.sha256:
                raise SystemExit(f"{p} changed since download (sha256 mismatch)")
            rows.append(dict(dataset=r.source, path=str(p.relative_to(ROOT)), size_bytes=size, sha256=sha, md5=md5,
                             verification=f"{r.verification}; url={r.url}; retrieved {r.retrieved}"))
    for p in sorted((ROOT / "data" / "external" / "nextgen_2026").glob("*")):
        sha, md5, size = hashes(p)
        rows.append(dict(dataset="nextgen_2026_figshare_29472362", path=str(p.relative_to(ROOT)), size_bytes=size,
                         sha256=sha, md5=md5, verification="md5 asserted in identify_release.py"))
    OUT.parent.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
    print(f"wrote {OUT} ({len(rows)} files)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
