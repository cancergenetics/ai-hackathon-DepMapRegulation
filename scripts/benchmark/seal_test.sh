#!/usr/bin/env bash
# Seal the held-out TEST labels (DESIGN v3.1 §6). RUN BY THE PI ONLY, once, after the freeze.
#
# 1. records sha256 of the frozen files (split, benchmark.py, evaluation scripts, universe/flag/annotation
#    tables, DepMap expression matrix, TRAIN labels, variants.yaml) and of every file to be sealed in results/FROZEN.sha256;
# 2. moves out of the project tree, to $DEPMAPTR_SEALED (default ~/DepMapTR_sealed), keeping relative paths:
#      data/processed/bench/positives_*_TEST.parquet, regulator_labels_*_TEST.tsv   (TEST labels)
#      data/external/pie_x_atlas_orion/de/*.parquet                                 (raw PIE DE tables)
#      every cache in data/processed/pie/ except the label-free allowlist below
#      (tested_*.parquet, loose_*.parquet, targets_*.tsv, target_summary_*.tsv hold per-perturbation rows,
#       i.e. TEST-regulator rows);
# 3. removes all read/write permissions on the sealed files (chmod a-rw); the sealed root is chmod 700
#    (sub-directories keep their mode). This is an audit-trail speed bump, not a security barrier: the
#    owner can chmod back, and the PIE DE tables are public. FROZEN.sha256, LOCK and the tuning log are
#    the actual safeguard.
#
# Allowlist kept in place (no per-perturbation rows, no labels): vocab_targets.tsv, vocab_features.tsv
# (identifiers only), features_<line>.tsv (per measured gene: control mean, number of perturbations tested).
# scripts/final_eval.py restores read access to the TEST label files for the single final run.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SEALED="${DEPMAPTR_SEALED:-$HOME/DepMapTR_sealed}"
FROZEN="$ROOT/results/FROZEN.sha256"
cd "$ROOT"

die() { echo "seal_test.sh: $*" >&2; exit 1; }

[ -f results/split/regulator_split.tsv ] || die "results/split/regulator_split.tsv missing: build the benchmark first"
[ -f data/processed/bench/manifest.json ] || die "benchmark tables missing"
if [ -e "$SEALED" ] && [ -n "$(ls -A "$SEALED" 2>/dev/null || true)" ]; then die "$SEALED is not empty: already sealed?"; fi
if [ -f "$FROZEN" ] && grep -q '^# sealed' "$FROZEN"; then die "$FROZEN already records a seal"; fi

shopt -s nullglob
TO_SEAL=(data/processed/bench/positives_*_TEST.parquet data/processed/bench/regulator_labels_*_TEST.tsv
         data/external/pie_x_atlas_orion/de/*.parquet)
for f in data/processed/pie/*; do
  case "$(basename "$f")" in
    vocab_targets.tsv|vocab_features.tsv|features_*.tsv) ;;
    *) TO_SEAL+=("$f") ;;
  esac
done
[ ${#TO_SEAL[@]} -ge 4 ] || die "expected TEST label files and raw DE tables; found: ${TO_SEAL[*]:-none}"
for f in "${TO_SEAL[@]}"; do [ -f "$f" ] || die "not a regular file: $f"; done

FROZEN_FILES=(results/split/regulator_split.tsv results/split/train_folds.tsv results/split/arm_units.tsv
              src/depmaptr/benchmark.py scripts/benchmark/build_benchmark.py
              scripts/final_eval.py scripts/evaluate.py scripts/benchmark/seal_test.sh
              data/processed/frozen_flags.tsv data/processed/panel_s2.tsv
              data/processed/gene_annotation.tsv data/processed/pie_gene_map.tsv
              data/processed/depmap_expr_matched.parquet
              data/processed/bench/universe_regulators_*.tsv data/processed/bench/universe_targets_*.tsv
              data/processed/bench/positives_*_TRAIN.parquet data/processed/bench/manifest.json
              results/variants.yaml)

{
  if [ -f "$FROZEN" ]; then cat "$FROZEN"; fi   # keep earlier entries (e.g. the split hash recorded at freeze)
  echo "# FROZEN.sha256 — block written by scripts/benchmark/seal_test.sh on $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "# git $(git rev-parse HEAD 2>/dev/null || echo nogit)"
  echo "# frozen (in tree; tests/test_frozen.py fails on drift)"
  for f in "${FROZEN_FILES[@]}"; do
    if [ -f "$f" ]; then shasum -a 256 "$f"; fi
  done
  echo "# sealed $(date -u +%Y-%m-%dT%H:%M:%SZ) -> $SEALED (path relative to the sealed root)"
  for f in "${TO_SEAL[@]}"; do shasum -a 256 "$f"; done
} > "$FROZEN.tmp"
[ -f results/variants.yaml ] || echo "seal_test.sh: WARNING results/variants.yaml not found (not hashed)" >&2

mkdir -p "$SEALED"
chmod 700 "$SEALED"
for f in "${TO_SEAL[@]}"; do
  mkdir -p "$SEALED/$(dirname "$f")"
  mv "$f" "$SEALED/$f"
  chmod a-rw "$SEALED/$f"
  [ ! -e "$f" ] || die "failed to move $f"
done
mv "$FROZEN.tmp" "$FROZEN"

# verify: nothing readable is left behind
for f in "${TO_SEAL[@]}"; do
  [ ! -e "$f" ] || die "$f still in the tree"
  if cat "$SEALED/$f" >/dev/null 2>&1; then die "$SEALED/$f is still readable"; fi
done
echo "sealed ${#TO_SEAL[@]} files into $SEALED; hashes in results/FROZEN.sha256 (commit it)."
