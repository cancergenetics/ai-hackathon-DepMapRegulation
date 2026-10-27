# Interfaces between pipeline components (PI-defined, 2026-10-27)

These contracts let the model engine (`assoc.py` and friends) and the benchmark harness (`benchmark.py`)
be built in parallel. Changes need a PI decision note.

## Gene identifiers
- The canonical gene key is the DepMap column label `"SYMBOL (ENTREZ)"`, e.g. `"TP53 (7157)"`.
- Chromosome, arm, start/end, PAR flag, and the PIE symbol mapping come from
  `data/processed/gene_annotation.tsv` and `data/processed/pie_gene_map.tsv` (data agent).

## Sample panel and frozen flags (data task T3 → data agent)
- `data/processed/panel_s2.tsv`: one row per S2 model, with ModelID and all covariates in raw form:
  - lineage, primary disease, sex;
  - ploidy, WGD (+ missing flags);
  - MSIScore (+ missing flag), TP53 status (+ missing flag);
  - media base, serum-free flag, growth pattern;
  - Stranded, platform;
  - PatientID, in_S1 flag.
- `data/processed/frozen_flags.tsv`: one row per gene key, with:
  - chrom, arm;
  - `variable_cn` (bool; X genes use the sex-centred test);
  - `measured` (bool);
  - `has_expr`, `has_cn`;
  - `is_par`, `is_y`, `is_mt`.

## Score files (model engine → harness)
`data/results/scores/<run_id>/` contains:
- `scores.npz`:
  - `regulators`: str[R], gene keys. **All** frozen variable-CN regulators (autosomes + X), not only
    perturbed genes. The tool does not depend on the benchmark.
  - `targets`: str[T], gene keys of all `measured` genes (excluding Y/PAR/MT).
  - `t`: float32[R,T], the signed reduced-form statistic (NaN allowed: no score).
  - `score`: float32[R,T], the signed final score S_AB (DESIGN §5.6). Ranking uses |score|; NaN → minimum.
  - `s_cis`: int8[R], s_A ∈ {−1, +1} (DESIGN §1).
  - `cis_t`: float32[R].
  - Same-chromosome entries may hold any value; the harness masks them.
  - **Optional fields** (D6; the harness uses them if present):
    - `call`: bool[R,T], the tool's call set (BH(z′) ≤ 0.05, or the MVP rule), used for the call-set OR.
      If it is absent, the OR is skipped with a note.
    - `pair_score`: float32[R,T], the signed *pre-fusion* pair component of `score` (I1: sign(t)·ρ_P,
      including any per-regulator factor such as I2's ρ_R).
    - `target_weight`: float32[T], the target-only fusion factor (I1: ρ_T^γ).

    These two must satisfy `score == pair_score * target_weight[None, :]`. If both are present, the C1
    control (depmaptr.controls) shuffles `pair_score` within bins and recomputes the fused score. If
    either is missing, it shuffles `score` directly, which is an approximation when `target_weight` varies
    within a C1 bin.
- `meta.json`: run_id, variant id, params, git hash, input checksums, n_samples, wall time, threads,
  peak RSS.

The predicted knockdown sign used by the harness is `−sign(score) · s_cis`.

## Benchmark tables (harness)
- `data/processed/bench/universe_regulators_{line}.tsv`: one row per regulator in that line's universe:
  - gene key, arm, split (TRAIN/TEST);
  - fold (1/2, TRAIN only);
  - own-gene class (i/ii/iii);
  - CN-variability tercile;
  - benchmark-degree decile (evaluation-only).
- `data/processed/bench/universe_targets_{line}.tsv`: one row per target in the line's universe, with
  gene key, chrom, 5 Mb bin and DepMap expression decile.
- `data/processed/bench/positives_{line}_{split}.parquet`: sparse rows for benchmark-tested trans pairs
  inside the universe with p_adj < 0.05, or meeting any sensitivity criterion. Columns: regulator,
  target, p_adj, lfc_sign, ueffect, and the flags pos_primary, pos_ueffect, pos_p01.
- The universe of a line is {regulators} × {targets} minus same-chromosome pairs and minus A = B.
  Pairs not in the positives file are negatives.
- TEST files are sealed after freeze (DESIGN §6). Tuning code must never open them.
