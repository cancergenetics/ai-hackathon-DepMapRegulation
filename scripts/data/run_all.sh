#!/usr/bin/env bash
# Rebuild the data foundation end-to-end (about 3 minutes plus downloads; peak RAM about 12 GB in cache_pie_de.py).
# Prerequisites: starting_data/ (DepMap 26Q3 files + HGNC protein-coding list) and the PIE DE tables in
# data/external/pie_x_atlas_orion/ (HF arcinstitute/PIE_x_atlas_orion @ abfa3da0f3bb98c9df9c90082386ed0c7c9282d3).
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=".venv/bin/python -I"
$PY scripts/data/download_depmap_extra.py      # 24Q4 figshare covariates, Sanger CMP, verify manual 26Q3 files
$PY scripts/data/download_annotation.py        # HGNC complete set, Ensembl GTF, UCSC cytobands
$PY scripts/data/identify_release.py           # release evidence -> data/processed/release_evidence.json
$PY scripts/data/convert_depmap.py             # CSV -> parquet (data/interim)
$PY scripts/data/cache_pie_de.py               # compact PIE caches (data/processed/pie)
$PY scripts/data/build_gene_annotation.py      # gene_annotation.tsv, pie_gene_map.tsv
$PY scripts/data/build_depmap_tables.py        # depmap_models.tsv, matched matrices, derived arm CN
$PY scripts/data/profile_depmap.py > /dev/null # depmap_profile.json
$PY scripts/data/profile_benchmark.py > /dev/null  # benchmark_profile.json, target summaries
$PY scripts/data/checksum_inputs.py            # docs/data_checksums.tsv
.venv/bin/python -m pytest -q tests/test_data_foundation.py
