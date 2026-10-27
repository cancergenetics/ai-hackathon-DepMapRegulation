#!/usr/bin/env python
"""Download covariate files that are not in starting_data/ and record URLs + checksums.

Sources (see DATA.md for the rationale):

1. DepMap Public 26Q3 (same release as starting_data/) -- the DepMap portal download API
   (https://depmap.org/portal/api/download/files) sits behind a human-verification page, so it
   cannot be scripted. Files from 26Q3 must be downloaded by a human in a browser from
   https://depmap.org/portal/data_page/?tab=allData and placed in data/external/depmap_26Q3/.
   This script then only *verifies and records* them (sha256/md5/size). Wanted files:
   OmicsGlobalSignatures.csv (ploidy, WGD, MSIScore, CIN, LoHFraction, Aneuploidy, PAM50) and,
   if present, OmicsArmLevelCNA.csv.

2. DepMap Public 24Q4 (figshare+ article 27993248, the last release with a figshare deposit and the
   last release in which PureCN ploidy / WGD / LoH / absolute CN were produced for new samples; per
   the 25Q2 release notes these features exist "only for samples released till 2024 Q4").
   Downloaded via the public figshare API, md5 verified against figshare's supplied_md5.
   Also the 24Q4 somatic mutation calls (MAF-like table + damaging/hotspot matrices) for TP53 status:
   26Q3 mutation files are portal-only, so 24Q4 is the latest programmatically available release.

3. Sanger Cell Model Passports model list (independent MSI status + ploidy annotation), mapped to
   DepMap via Model.csv SangerModelID.

Usage:  .venv/bin/python scripts/data/download_depmap_extra.py [--skip-absolute]
Writes: data/external/<source>/..., data/external/download_manifest.tsv
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import truststore  # use the OS certificate store (python.org builds lack CA certs)

truststore.inject_into_ssl()

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "data" / "external"
MANIFEST = EXT / "download_manifest.tsv"

FIGSHARE_24Q4_ARTICLE = 27993248
FILES_24Q4 = [
    "README.txt",
    "OmicsSignatures.csv",
    "OmicsSignaturesProfile.csv",
    "OmicsDefaultModelProfiles.csv",
    "OmicsProfiles.csv",
    "Model.csv",
    "OmicsAbsoluteCNSegmentsProfile.csv",
    "OmicsAbsoluteCNGene.csv",
    # somatic mutation calls (Mutect2) for TP53 functional status (DESIGN v2 data task T1)
    "OmicsSomaticMutationsMatrixDamaging.csv",
    "OmicsSomaticMutationsMatrixHotspot.csv",
    "OmicsSomaticMutations.csv",
]
ABSOLUTE_FILES = {"OmicsAbsoluteCNGene.csv"}  # large (239 MB); optional

CMP_URL = "https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz"

MANUAL_26Q3 = ["OmicsGlobalSignatures.csv", "OmicsArmLevelCNA.csv"]


def hashes(path: Path) -> tuple[str, str, int]:
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            sha.update(chunk)
            md5.update(chunk)
    return sha.hexdigest(), md5.hexdigest(), path.stat().st_size


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as out:
        while True:
            chunk = r.read(1 << 22)
            if not chunk:
                break
            out.write(chunk)
    tmp.rename(dest)


def figshare_files(article: int) -> dict[str, dict]:
    url = f"https://api.figshare.com/v2/articles/{article}/files?page_size=500"
    with urllib.request.urlopen(url, timeout=120) as r:
        return {f["name"]: f for f in json.load(r)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-absolute", action="store_true", help="skip OmicsAbsoluteCNGene.csv (239 MB)")
    args = ap.parse_args()
    rows = []
    today = dt.date.today().isoformat()

    # --- 2. DepMap 24Q4 via figshare -------------------------------------------------------
    d24 = EXT / "depmap_24Q4"
    listing = figshare_files(FIGSHARE_24Q4_ARTICLE)
    for name in FILES_24Q4:
        if args.skip_absolute and name in ABSOLUTE_FILES:
            continue
        meta = listing[name]
        dest = d24 / name
        if not dest.exists() or dest.stat().st_size != meta["size"]:
            print(f"downloading {name} ({meta['size']/1e6:.1f} MB)", file=sys.stderr)
            fetch(meta["download_url"], dest)
        sha, md5, size = hashes(dest)
        if md5 != meta["supplied_md5"] or size != meta["size"]:
            raise SystemExit(f"checksum mismatch for {name}: {md5} vs {meta['supplied_md5']}")
        rows.append(("depmap_24Q4", name, meta["download_url"], f"figshare:{FIGSHARE_24Q4_ARTICLE}",
                     size, sha, md5, "md5 matches figshare supplied_md5", today))

    # --- 3. Sanger Cell Model Passports ----------------------------------------------------
    dest = EXT / "sanger_cmp" / "model_list_latest.csv.gz"
    if not dest.exists():
        print("downloading Cell Model Passports model list", file=sys.stderr)
        fetch(CMP_URL, dest)
    sha, md5, size = hashes(dest)
    rows.append(("sanger_cmp", dest.name, CMP_URL, "model_list_latest (gzip header mtime records version)",
                 size, sha, md5, "no upstream checksum published", today))

    # --- 1. DepMap 26Q3 manual files -------------------------------------------------------
    d26 = EXT / "depmap_26Q3"
    d26.mkdir(parents=True, exist_ok=True)
    for name in MANUAL_26Q3:
        p = d26 / name
        if p.exists():
            sha, md5, size = hashes(p)
            rows.append(("depmap_26Q3", name, "https://depmap.org/portal/data_page/?tab=allData (manual)",
                         "DepMap Public 26Q3", size, sha, md5, "manual browser download; compare md5 to portal listing", today))
        else:
            print(f"NOTE: {p} absent -- portal requires browser verification; download manually if needed.",
                  file=sys.stderr)

    with open(MANIFEST, "w") as fh:
        fh.write("source\tfile\turl\tversion\tsize_bytes\tsha256\tmd5\tverification\tretrieved\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {MANIFEST} ({len(rows)} files)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
