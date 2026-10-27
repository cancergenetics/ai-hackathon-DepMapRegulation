#!/usr/bin/env bash
# Run a tuning tier: build each variant's scores and compare to the incumbent on HCT116 TRAIN (DESIGN §5.5, §8.2).
# usage: scripts/tune_tier.sh <tier-name> <incumbent-run> <variant-spec>...
#   variant-spec is passed to run_model.py --variant (e.g. "g1_t30_w0" or "focal+g1_t30_w0"); run id = spec with '+' -> '__'
set -euo pipefail
cd "$(dirname "$0")/.."
tier=$1; inc=$2; shift 2
mkdir -p results/tuning
runs=()
for v in "$@"; do
  rid=${v//+/__}
  if [ ! -f data/results/scores/$rid/scores.npz ]; then
    .venv/bin/python scripts/run_model.py --variant "$v" --out data/results/scores/$rid > results/tuning/run_$rid.log 2>&1
  fi
  runs+=("$rid")
done
.venv/bin/python scripts/evaluate.py "${runs[@]}" --compare "$inc" --line hct116 > results/tuning/${tier}_compare.json 2>&1
.venv/bin/python scripts/summarise_compare.py results/tuning/${tier}_compare.json
