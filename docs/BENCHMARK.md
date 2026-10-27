# Benchmark harness (DESIGN v3.1 §6–§8)

Code: `src/depmaptr/benchmark.py` (library), `scripts/benchmark/build_benchmark.py` (tables and split),
`scripts/evaluate.py` (TRAIN evaluation CLI), `scripts/final_eval.py` (the single TEST run, plus `--mock`),
`scripts/benchmark/seal_test.sh` (run by the PI only), `scripts/benchmark/calibrate_b0.py`, and
`scripts/benchmark/split_cn_diagnostic.py` (label-free split-leakage diagnostic).
Tests: `tests/test_benchmark.py` (toy known answers), `tests/test_final_eval.py` (seal → final_eval
workflow on a throwaway toy git repo), `tests/test_benchmark_real.py` (real-data invariants), and
`tests/test_frozen.py` (hash drift, active once `results/FROZEN.sha256` exists).

The benchmark labels are used only inside the evaluation functions. Nothing in this module produces or
modifies a score.

## 1. Tables (built 2026-10-27 from frozen_flags.tsv and panel_s2.tsv by the data agent)

| | HCT116 (primary) | HEK293T |
|---|---:|---:|
| PIE perturbations (`target`) | 8,987 | 4,297 |
| − unmapped to a DepMap key (incl. 38 symbol collisions, map-wide) | 61 | 56 |
| − Y / PAR / MT / unplaced | 16 | 12 |
| − not frozen variable-CN | 3,287 | 1,541 |
| **universe regulators** | **5,623** (360 on X; 3 CN-only) | **2,688** (180 on X; 3 CN-only) |
| consistently tested features (control mean ≥ 0.05) | 11,100 | 12,247 |
| − unmapped | 728 | 1,072 |
| − Y / PAR / MT | 24 | 24 |
| − not frozen `measured` | 139 | 151 |
| **universe targets** | **10,209** (572 × 5 Mb bins) | **11,000** (584 bins) |
| universe pairs (trans, A ≠ B) | 54,357,451 | 27,984,730 |
| primary positives (p_adj < 0.05) | 115,601 | 52,510 |
| median on-target LFC (own feature tested; TRAIN regulators) | −27.6 | −26.7 |

Every consistently tested feature is tested for every perturbation (asserted at build time), so every
universe pair is a tested pair.

## 2. Held-out split (arm level; seed 20261027)

There are 42 arm units: all arms of chr1–22 and X that hold ≥ 1 universe regulator in either line,
including 2p (2 HCT116 regulators) and 21p (1).

**TEST arms (14):** 17q, 7q, 14q, 6q, 3p, Xp, 7p, 22q, 19p, 4q, 5p, 12q, 18p, 21q.

The 28 TRAIN arms alternate folds 1/2 in sorted order (`results/split/train_folds.tsv`).
`results/split/arm_units.tsv` lists every unit with its order, triplet, split and fold.
`results/split/regulator_split.tsv` gives the gene → arm → split → fold mapping, with line membership.

Positive totals for TEST are aggregate split sizes. Other label-derived TEST summaries (own class,
direction, per-regulator degree) are deliberately not reported anywhere in the tree.

| line | split | arms | regulators | pairs | K at 0.1 % | pos. p_adj<0.05 | pos. \|Ueffect\| ≥ 0.2 (D4) | pos. p<0.01 | % LFC<0 | own class i / ii / iii |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| HCT116 | TRAIN | 28 | 3,806 | 36,664,698 | 36,665 | 81,258 | 38,342 | 42,581 | 58.9 | 1,388 / 171 / 2,247 |
| HCT116 | TEST | 14 | 1,817 (32 %) | 17,692,753 | 17,693 | 34,343 | 14,748 | 18,391 | — | — |
| HEK293T | TRAIN | 28 | 1,890 | 19,611,225 | 19,611 | 40,343 | 6,232 | 14,131 | 44.5 | 147 / 38 / 1,705 |
| HEK293T | TEST | 14 | 798 (30 %) | 8,373,505 | 8,374 | 12,167 | 3,586 | 3,460 | — | — |

TRAIN folds:
- HCT116: fold 1 has 14 arms, 1,989 regulators and 42,454 positives; fold 2 has 14 arms, 1,817 regulators and 38,804 positives.
- HEK293T: fold 1 has 1,006 regulators and 20,881 positives; fold 2 has 884 regulators and 19,462 positives.

2,039 regulators are in both lines' universes.

## 3. Exact definitions

- **Gene key:** the DepMap column label. A PIE id maps to the `depmap_cn_column`, else to the
  `depmap_expr_column` of `pie_gene_map.tsv`. The mapping is asserted 1:1.
- **Regulators of a line:** PIE `target` in that line's DE table, mapped, `variable_cn` (frozen flags),
  chromosome 1–22 or X, not PAR/Y/MT.
- **Targets of a line:** PIE `feature` with control mean ≥ 0.05 in that line (`features_<line>.tsv`
  `ref_mean`), mapped, `measured`, chromosome 1–22 or X, not PAR/Y/MT.
- **Universe:** regulators × targets minus same-chromosome pairs, which also removes A = B.
- **Positives:** universe pairs with PIE `p_adj < 0.05`, built from `loose_<line>.parquet`, which holds
  every row with p_adj < 0.1 at full precision. Each positive carries these flags:
  - `pos_primary`: p_adj < 0.05;
  - `pos_ueffect`: additionally |Ueffect| ≥ 0.2 (**PI decision D4**). PIE `Ueffect` is a signed
    rank-biserial effect in [−1, 1], so the v3.1 text |Ueffect − 0.5| ≥ 0.1 (written for a [0, 1] scale)
    kept 98.6 % of positives. |Ueffect| ≥ 0.2 is the equivalent of |AUC − 0.5| ≥ 0.1. It keeps 47 %
    (HCT116 TRAIN) and 15 % (HEK293T TRAIN) of the positives;
  - `pos_p01`: p_adj < 0.01;
  - `lfc_sign` = sign(LFC). Only the sign is stored; the magnitude is never used.
- **Own-gene class** (from the PIE row whose feature is the perturbed gene itself, same PIE id):
  - (i) tested and p_adj < 0.05 and LFC < 0;
  - (ii) tested otherwise;
  - (iii) not tested.

  Tested status comes from `tested_<line>.parquet`, and p_adj from `loose`. A tested own row absent from
  `loose` has p_adj ≥ 0.1, so it is class (ii).
- **CN-variability tercile:** terciles, among the line's universe regulators, of the frozen flag statistic
  `n_lines_abs_log2cn_ge1`. This is the number of S2 lines with |log2 relCN| ≥ 1, sex-centred for X: the
  statistic that defines `variable_cn`. HCT116 cut points are 20–28, 29–50 and 51–785.
- **Benchmark-degree decile** (evaluation only): the number of primary positives per regulator, in deciles
  computed **within split**.
- **Quantile bins** (terciles and deciles): bin = ⌊n·(average rank − ½)/N⌋ + 1. Tied values always share a
  bin. About 52 % of TRAIN regulators have degree 0, so the degree "deciles" collapse to 5 occupied levels
  (3, 7, 8, 9, 10). This is harmless for stratification.
- **Target strata:**
  - the 5 Mb bin is `chrom:⌊gene midpoint / 5e6⌋`;
  - the expression decile is the decile of the mean log2(TPM+1) over panel S2, among the line's universe
    targets.
  - The table also carries DepMap expression in HCT116 (ACH-000971) and HEK293 (ACH-001085), for B4e/NC2.
  - Regulators carry `expr_context`, their DepMap expression in the line's context model (label-free; for
    G2/NC2).
- **Split:**
  - The arm units are sorted by the number of HCT116 universe regulators, **descending**, with ties in
    genomic order.
  - Consecutive triplets are drawn with `np.random.default_rng(20261027).integers(3)`, one draw per
    triplet in order. The drawn position is TEST.
  - A final partial group (none here: 42 = 14 × 3) would get a TEST arm only if the draw is smaller than
    its size.
- **Score handling:**
  - Universe genes absent from a score file get NaN, which counts as the minimum score.
  - Ranking uses |score|, with NaN at the minimum.
  - Same-chromosome entries are masked whatever value they hold.
  - Ties at the top-K boundary are broken by a seeded hash: blake2b of `seed|gene key` per gene, combined
    by splitmix64. The result does not depend on row order.
  - K = round(frac × universe pairs of the evaluated split and line).

## 4. Metrics (all per line × split; `Evaluation` / `compute_metrics`)

**Pooled overlap**
- Pooled 2×2 at K ∈ {0.1 %, 0.01 %, 1 %}: sample OR ad/bc, plus the Fisher p (one-sided "greater" and
  two-sided). The Fisher p is labelled "assumes independence".
- When a cell is 0, every OR in this module uses the Haldane +0.5 correction.

**Regulator-label permutation** (1,000 permutations)
- Label rows are reassigned among the split's regulators within strata of own class × CN tercile.
- Replicate k uses the pairs trans to both A and π(A). The observed OR of that replicate uses the same
  pairs, and the top-K selection stays fixed.
- p = (1 + #{OR_null ≥ OR_obs}) / (n + 1).
- Also reported: the median null OR and OR / median null OR.
- The secondary permutation uses strata of own class × degree decile.
- The computation is O(R) per replicate, using M = Sel·Posᵀ (R×R) and per-regulator × target-chromosome
  counts. It matches a dense brute force exactly (test).

**Jackknife and bootstrap**
- Delete-one-arm jackknife of log OR, with targets fixed, top-K fixed and t_{G−1} critical values. G = 28
  (TRAIN) or 14 (TEST).
- The two-way bootstrap is a sensitivity analysis: multinomial weights on regulator arms × target 5 Mb
  bins, applied to the per-(arm × bin) 2×2 counts. It reports a percentile CI.

**CMH**
- Mantel–Haenszel OR with strata regulator × target-expression decile, and own class × expression decile.
- Also reported: the Robins–Breslow–Greenland SE and the CMH χ² p (no continuity correction).
- Both match statsmodels `StratifiedTable`.

**Other secondary metrics**
- Per-regulator AUROC: the Mann–Whitney AUC of |score| with average ranks (NaN at the minimum) over the
  regulator's trans targets. The mean is over regulators with ≥ 10 positives.
- Share of top-K pairs from the 10 regulators with the most top-K pairs.
- Hub-excluded sensitivity:
  - drop the top 1 % of regulators and of targets by degree within the evaluated split (hash tie-break);
  - **re-rank top-K within the reduced universe**.
- Results for the own-class (i) and (iii) sub-universes (top-K re-ranked within each).
- Call-set OR, from an optional boolean `call` array [R, T] in `scores.npz` (aligned with `regulators` ×
  `targets`; the model side will add it) or the `call_mask` argument. Without one, `metrics.callset` holds
  `{"skipped": ...}`. Note: in the sensitivity-set ORs, primary positives that fail the stricter criterion
  count as negatives (they are not excluded from the universe).

**Sign metrics** (primary positives; top 0.1 %)
- Predicted knockdown sign = −sign(score) · s_cis. It is 0 when the score is NaN or 0, and 0 never counts
  as concordant.
- Sign-concordant OR: positives are the pairs that are benchmark-positive and concordant.
- Among the overlaps (top-K ∩ positives):
  - the concordance rate;
  - c₀ = q·u + (1−q)(1−u), where q is the fraction predicted down and u the fraction with LFC < 0;
  - the 2×2 sign OR (pred down/up × LFC down/up);
  - the within-regulator sign permutation p (LFC signs shuffled within each regulator's overlaps);
  - regulator-weighted concordance (regulators with ≥ 5 overlaps).
- A unit test covers a regulator with s_cis = −1.

**Controls built in the harness** (evaluation-only references; never used for tuning):
- **PC2** (`evaluate_pc2`, `evaluate.py --pc2`; also run by `final_eval.py`):
  - universe: HCT116 regulators of the split that are also HEK293T universe regulators (TRAIN: 1,441; TEST: 598), ×
    all HCT116 universe targets;
  - score of (A, B) = (HEK293T LFC sign) × −log10(HEK293T p_adj) for HEK293T primary positives, else 0;
  - s_cis = +1, so the predicted knockdown sign is the HEK293T LFC sign;
  - `call` = HEK293T positive, so `callset.or` is the plain overlap OR of the two hit sets. The top-K ORs
    rank HEK293T hits by significance; HCT116 targets not in the HEK293T universe score 0.
- **NC2** (subset `nc2`; `evaluate.py --subset nc2`; a `subsets` entry in `selected.yaml`) is the
  regulators with own-gene class (iii) and DepMap context expression log2(TPM+1) < 1.
  - Context: on the HCT116 line, DepMap HCT116 (ACH-000971), as DESIGN §7 states. On the HEK293T line,
    HEK293 (ACH-001085), which is the analogous context (DESIGN names only HCT116).
  - The 3 CN-only regulators per line (no expression column) count as not expressed, as in DESIGN §4.
  - NC2 must be run on a score file **without G2**. The harness cannot check this; `selected.yaml` must
    point to the no-G2 run.
  - Expected OR ≈ 1.
- **Split-sibling sensitivity** (subset `no_split_sibling`; run automatically by `final_eval.py` for the
  selected method) keeps the regulators whose arm's sibling arm is not in the other split (§9).

**Runtime.** `evaluate()` logs its own wall time, and copies the score run's wall time and threads from
`meta.json`. A full HCT116 TRAIN evaluation (1,000 permutations and 1,000 bootstraps) takes about 4 s.

**Paired comparison** (`compare()`, `evaluate.py --compare`):
- Δlog OR per TRAIN fold;
- the paired arm jackknife;
- the paired two-way bootstrap, with shared weights: the Δlog OR 90 % lower bound and the Δlog CMH OR
  90 % upper bound;
- `adopt_rule_8_2` = both fold gains > 0, the LB > 0 and the CMH UB > 0.

**Fold ORs.** The fold ORs use the full-TRAIN top-K selection restricted to each fold's regulators. This
is the same "threshold fixed at its full-split value" convention as the resampling replicates.

**Seeds.** Seeds come from `numpy.random.SeedSequence(20261027).spawn(8)` at these spawn indices:
- 0: regulator permutation;
- 1: degree permutation;
- 2: sign permutation;
- 3: bootstrap.

The B0 calibration uses spawn indices 100+. The tie-break hash uses seed 20261027, and the split uses
`default_rng(20261027)`.

## 5. Guards, tuning log, sealing

**`evaluate()` and `compare()`**
- They raise `SealedSplitError` for any split other than TRAIN, unless they are given the token that only
  `final_eval_token()` returns. That function requires the running program to be `scripts/final_eval.py`.
- `read_positives()` also refuses any file whose name contains TEST.

**Tuning log.** Every `evaluate()` and `compare()` call appends a row. Logging cannot be switched off;
the old `log=False` was removed. `debug=True` (`evaluate.py --debug`; smoke tests, harness checks) sends
the row to `results/tuning_log_debug.tsv` instead, which is also committed. Each row holds:
- run id, variant, kind, line and split;
- the params hash (meta params plus evaluation parameters);
- the git hash with a dirty flag;
- data checksums (sha256 of the bench tables used and of `scores.npz`);
- the seed and spawn map, n_perm and n_boot;
- the headline metrics as JSON;
- the wall time.

Compare rows carry both variants and both params hashes. Full metrics go to
`results/eval/<run>/<line>_<split>.json` (debug runs: `results/eval/debug_<run>/`).

Remaining bypasses (by design, a speed bump):
- calling `Evaluation` methods directly, as `calibrate_b0.py` and the real-data tests do with random
  scores;
- pointing `DEPMAPTR_ROOT` elsewhere.

**The audit trail is the actual safeguard:** `FROZEN.sha256`, the LOCK commit and the tuning logs.

**TEST-only label columns.** `own_class` and `degree_decile` of TEST regulators are not in the open
`universe_regulators_<line>.tsv`, where they are empty. They are in `regulator_labels_<line>_TEST.tsv`,
which is sealed with the TEST positives (implementation choice; see open issue 4). `manifest.json` and
`results/split/benchmark_summary.json` carry label-derived summaries (own class, median on-target LFC,
direction, regulators with hits) for TRAIN only. For TEST they carry only arms, regulators, pairs and
positive totals.

**`seal_test.sh`** (the PI runs it; it has only been exercised on a toy copy in the tests):
1. It appends the sha256 of the frozen files and of every file it seals to `results/FROZEN.sha256`. The
   frozen files are:
   - the split files;
   - `benchmark.py`, `build_benchmark.py`, `final_eval.py`, `evaluate.py` and `seal_test.sh`;
   - the flags, panel, `gene_annotation.tsv`, `pie_gene_map.tsv` and `depmap_expr_matched.parquet`;
   - the universe tables, the TRAIN positives, the manifest and `variants.yaml`.
2. It moves these to `~/DepMapTR_sealed/<same relative path>` and runs `chmod a-rw` on them:
   - the TEST positives and regulator-label files;
   - `data/external/pie_x_atlas_orion/de/*.parquet`;
   - everything in `data/processed/pie/` except the label-free `vocab_*.tsv` and `features_<line>.tsv`.
     That is: `tested_*`, `loose_*`, `targets_*`, `target_summary_*`.
3. It refuses to run twice.

After sealing, `build_benchmark.py` cannot be rerun, by design.

**`final_eval.py`** (real run)
1. Guards. It refuses if:
   - `results/final/LOCK` exists in the tree or in HEAD;
   - `results/final/selected.yaml` is missing, untracked or modified;
   - `results/FROZEN.sha256` is missing, or any in-tree frozen file has drifted (it checks this itself;
     it does not rely on pytest);
   - TEST result JSON already exists under `results/final/`;
   - `tuning_log.tsv` already has TEST rows;
   - more than 3 ablations are listed, a subset is unknown, or a listed run directory is missing.
2. It writes `results/final/LOCK` (status STARTED, timestamp, `selected.yaml` sha256, git hash) and
   **commits it before any TEST file is read**. A crash after this point therefore still counts as the
   one run.
3. It sets read access (u+r) on the TEST label files only, and checks each one's sha256 against the
   sealed block of `FROZEN.sha256`. A mismatch fails the run.
4. It evaluates:
   - every listed run on TEST, for both lines;
   - the `subsets` (e.g. NC2);
   - `no_split_sibling` for the selected run;
   - PC2;
   - the paired comparisons of the selected run against each other run.
5. In a `finally` block, it removes read access. LOCK becomes COMPLETED, or FAILED with the error, and
   `results/final/` is committed again.

**`--mock`** (the DESIGN §6 dry run) goes through the **same** code path:
- `build_mock_sealed` writes TRAIN fold-2 positives and labels into `data/interim/mock_sealed/` under the
  real TEST file names, plus its own `FROZEN.sha256`, and sets them `chmod a-rw`;
- the mock token's `load_bench` takes the TEST branch: sealed-dir reads, the 1:1 label merge, hash
  checks, the chmod round trip and the LOCK lifecycle;
- output goes to `results/final_mock/` (split PSEUDO_TEST), with no git commit;
- it can be repeated.

The workflow is tested end to end on a toy git repo in `tests/test_final_eval.py`. The tests cover the
mock, the seal, drift refusal, refusal of an uncommitted `selected.yaml`, the LOCK commit order, and a
crash after the LOCK, which consumes the run.

## 6. Calibration on the real TRAIN universe under B0 (200 seeds, `calibrate_b0.py`)

| | HCT116 | HEK293T |
|---|---|---|
| regulator permutation p: mean / P(p ≤ .05) / KS p | 0.498 / 0.045 / 0.66 | 0.495 / 0.030 / 0.82 |
| degree-decile permutation p | 0.499 / 0.025 / 0.996 | 0.492 / 0.040 / 0.50 |
| Fisher p (one-sided) | 0.513 / 0.055 / 0.32 | 0.520 / 0.025 / 0.09 |
| CMH p | 0.494 / 0.055 / 0.81 | 0.506 / 0.025 / 0.74 |
| sign permutation p | 0.663 / 0.000 (conservative) | 0.691 / 0.005 (conservative) |
| jackknife 95 % CI coverage of log OR = 0 | 0.92 | 0.97 |
| two-way bootstrap 95 % CI coverage | 1.00 (conservative) | 1.00 (conservative) |
| median OR | 0.997 | 0.991 |

**Permutation tests.** The regulator-label permutation tests are calibrated, as DESIGN §6 requires.

**Toy-size limits.** On toys with only 4 chromosomes and tiny counts, the regulator permutation is
conservative (P(p ≤ .05) ≈ 0.01). Each replicate then drops a large share of the pairs, which adds noise
to the observed OR. The unit tests therefore use 23 chromosomes.

**Sign permutation.** Under B0 most regulators have about one overlap, so the within-regulator shuffle
cannot move. The test is valid but has little power in that regime.

**Jackknife coverage.** At 92 % for HCT116 (binomial SE about 2 %), it is slightly below nominal.

**NC1.** Calibration on NC1 is still needed (DESIGN §6). Run `calibrate_b0.py`-style loops on NC1 score
files once they exist.

## 7. Minor implementation choices (not in DESIGN; recorded here)

1. Arm sort is descending by HCT116 count, with ties in genomic order; the RNG draws one integer per
   triplet.
2. Quantile bins use average ranks, and ties share a bin.
3. The degree decile is computed within split.
4. The CN tercile uses the frozen `n_lines_abs_log2cn_ge1`.
5. The 5 Mb bin uses the gene midpoint. Seven DepMap genes have no coordinates, only an HGNC band (e.g.
   PAXX). They get their UCSC hg38 band extent.
6. If any 2×2 cell is 0, ORs use the Haldane correction.
7. K = round(frac·N).
8. AUROC uses average ranks.
9. Hub exclusion, the class (i)/(iii), NC2 and split-sibling sub-universes, and PC2 all re-rank the top-K
   within their universe. **Fold ORs keep the full-TRAIN top-K selection** restricted to each fold's
   regulators. That follows the §6 "threshold fixed at its full-split value" convention; DESIGN §8.2's
   "gain in both folds" would also allow a per-fold re-rank.
10. The 2×2 sign table excludes pairs with a predicted or LFC sign of 0.
11. Interface naming: DESIGN §6 mentions `bench_train.parquet` / `bench_test.parquet`, while
    INTERFACES.md specifies `positives_<line>_<split>.parquet`. The INTERFACES names are used.
12. The context-model assertion uses `StrippedCellLineName` (HCT116 / HEK293), because Model.csv spells
    `HEK-293`.
13. The NC2 context on the HEK293T line is HEK293; CN-only regulators count as not expressed (§4).

## 8. Open issues for the PI

1. **Ueffect.** Resolved by D4 (|Ueffect| ≥ 0.2; tables rebuilt by the PI and rebuilt again here with
   byte-identical outputs).
2. **NC1 calibration of the permutation test** has not been run (needs NC1 score files).
3. **Two-way bootstrap and jackknife.** The bootstrap is conservative (100 % coverage), and the
   jackknife is about 92 % on HCT116. This matters for the §8.2 adoption rule, which uses the paired
   bootstrap LB: the rule will be strict.
4. **TEST label columns.** `own_class` and `degree_decile` for TEST regulators are split into a sealed
   file instead of being kept in `universe_regulators_<line>.tsv` (INTERFACES lists them in one table).
5. **`reviews/round2/data.md`** names a TEST-arm hub regulator and its hit count. It was written before
   the split existed; it is not a harness file, so it was left unchanged for the PI to decide.

## 9. Split-arm leakage (DESIGN-conformant, a reported limitation)

The split unit is the arm, so for 10 chromosomes one arm is TEST and the other TRAIN: 3, 4, 5, 6, 12, 17,
18, 19, 21 and X. 7p/7q are both TEST; 14q, 22q and 21q's partners are absent or tiny. A whole-chromosome
CN event makes regulators on sibling arms share CN vectors, which is a weaker form of the leakage DESIGN
gave for rejecting 5 Mb bins.

`scripts/benchmark/split_cn_diagnostic.py` uses only DepMap CN and the split, no labels. It writes
`results/split/cn_leakage_*.tsv` and `cn_leakage_summary.json`.
- **Positive control:** within an arm, gene–gene CN correlation is high (8q genes: median r = 0.61);
  between 8q and 8p it is 0.04.
- **Regulators.** For each TEST regulator it records the max |r| with any TRAIN regulator:
  - HCT116: median 0.19 and 90th percentile 0.34. 7.4 % are > 0.5, and none are > 0.8.
  - These are larger for TEST regulators whose sibling arm is TRAIN: median 0.20, 90th percentile 0.58,
    12 % > 0.5. Without a TRAIN sibling: 0.16, 0.21 and 0 %.
  - For reference, TRAIN fold 2 vs fold 1 gives median 0.14 and 90th percentile 0.27.
  - HEK293T is similar.
- **Arms.** Arm-CN |r| between a TEST arm and its TRAIN sibling has median 0.19. The exception is
  **Xp–Xq, r = 0.55**, from whole-X loss (D3). Others: 4q–4p 0.29, 12q–12p 0.26, 18p–18q 0.22, 21q–21p
  0.20, 6q–6p 0.18, and ≤ 0.10 for 17, 3, 19 and 5.

The leakage is therefore modest except on X.
- The `no_split_sibling` subset gives an evaluation-only sensitivity analysis. It keeps the TEST arms
  whose sibling is TEST or absent: 7p, 7q, 14q and 22q (721 HCT116 TEST regulators). 21q is excluded because its sibling 21p is TRAIN, although 21p carries only
  1 regulator.
- REPORT.md should state the limitation and show this sensitivity analysis.

## 10. D6 matched controls and diagnostics (added after the freeze; no frozen metric changed)

Code: `src/depmaptr/controls.py`, `scripts/evaluate.py --c1 / --compare-c1 / --final-candidate`.
`scripts/final_eval.py` also runs C1 and reports the D6 CMH.

**C1 (pair shuffle)**
- **What is shuffled.** Each regulator's score row is permuted across targets within 50 bins of DepMap
  context-expression percentile:
  - the context is HCT116 (ACH-000971) on the HCT116 line, and HEK293 on the HEK293T line;
  - percentiles are taken over **all score-file targets** (all DepMap measured genes), not the PIE
    universe;
  - same-chromosome entries are permuted only among themselves.
- **Seed.** The order comes from seeded pair-hash keys: `SeedSequence(20261027).spawn(8)[4]` seeds
  blake2b/splitmix64. The result does not depend on row order or on which regulators are evaluated.
- **Exact mode.** If `scores.npz` has `pair_score` [R,T] and `target_weight` [T] (INTERFACES, optional;
  `score = pair_score · target_weight`), C1 shuffles `pair_score` and recomputes the fused score. This is
  exact for the I1/I2 fusion.
- **Approximate mode.** Otherwise C1 shuffles `score` directly. Results report `c1_mode` and carry a
  `c1_warning`.
  - **Known bias:** on HCT116 TRAIN (smoke test, debug log), B4e is a pure target-expression score, yet its
    OR goes from 7.09 to 4.08 under C1 (Δlog OR +0.55, jackknife CI 0.40 to 0.71). This happens because
    the top 2 % expression bin still has a steep detectability gradient, so the approximation overstates
    pair-specific gain for target-weighted scores.
  - **Fused runs must supply `pair_score` and `target_weight`.**
  - Smoke test of the incumbent (`g1_t100_w0__a3_r095_s2_w0`, approximate mode, which is exact for a
    score with no target factor): OR 5.00 vs 4.95 under C1, Δlog OR +0.01 (CI −0.22 to 0.24).
- **Pair-specific gain** (`compare_c1`, logged as kind `control:C1`):
  - OR(run) − OR(C1) and Δlog OR;
  - a paired delete-one-arm jackknife with t_{G−1};
  - a paired two-way bootstrap of Δlog OR (shared weights);
  - top-K is re-selected for each score.
- `evaluate_c1` runs the full metrics of the C1 score (`run_id` `<run>_C1`).

**D6 diagnostic CMH**
- Key `cmh_regulator_x_context_exprdecile_D6`, label "CMH_OR(reg x HCT116 expr) [D6 diagnostic]".
- Strata are regulator × deciles of DepMap context expression among the line's universe targets.
- The frozen regulator × S2-mean-expression CMH is unchanged.
- It reproduces the review's diag2 (B4e 3.38).

**Permutation counts.** `evaluate.py` already defaulted to `--n-perm 1000`; the 200s in the tree came
from explicit flags. `--final-candidate` enforces n_perm and n_boot ≥ 1000 and adds `--compare-c1`.

**Frozen-file change.** `benchmark.py`, `evaluate.py` and `final_eval.py` were edited additively. The new
`controls.py` is not yet frozen. These need a decision note, then:

    .venv/bin/python scripts/benchmark/record_frozen_change.py --note "D6 controls" \
        src/depmaptr/benchmark.py src/depmaptr/controls.py scripts/evaluate.py scripts/final_eval.py

Until then, `tests/test_frozen.py` fails, and `final_eval.py` refuses, as intended.

## 11. D8 evaluation-side sensitivity analyses (rules fixed before any result was computed)

### 11.1 Regulator classes (report analysis; HGNC gene groups)

**Source.** `starting_data/gene_with_protein_product.txt` (HGNC, PI-provided; checksum in DATA.md),
column `gene_group`. This holds the leaf group names, separated by `|`.

**Mapping.** DepMap gene key → `hgnc_id` (`gene_annotation.tsv`) → HGNC row. Each group name is tested
against the rules below; matching is case-sensitive on the exact HGNC names unless stated otherwise. A gene
in several classes takes the first one in the order **ribosomal > RBP > TF > kinase > other**. Genes with
no HGNC row or no group are "other".

| class | a gene group matches (regular expressions on single group names) |
|---|---|
| ribosomal | `ribosomal proteins$`: L and S ribosomal proteins, mitochondrial, large/small subunit mitochondrial. Not "Ribosomal protein S6 kinase family" or biogenesis complexes. |
| RBP | `RNA binding` (any group containing it, e.g. "RNA binding motif containing", ELAV, Musashi, DAZ, Mex-3, fox-1, FASTK); `^Heterogeneous nuclear ribonucleoproteins$`; `^DEAD-box helicases$`; `^DEAH-box helicases$` |
| TF | `^Zinc fingers C2H2-type$`; `homeobox` (case-insensitive: all homeobox classes and families); `^Basic helix-loop-helix proteins$`; `bHLH transcription factor`; `^Basic leucine zipper proteins$`; `bZIP transcription factor`; `^Forkhead boxes$`; `^Nuclear receptor subfamily`; `^GATA zinc finger domain containing$`; any group containing `transcription factor` (case-insensitive) except those starting `General transcription factor` |
| kinase | the group name contains `kinase` (case-insensitive substring, so it includes e.g. pseudokinase, hexokinase, lipid and metabolic kinases), unless it contains one of `anchoring`, `kinase phosphatase`, `kinase inhibitor`, `regulatory subunit`, `kinase activators`, `non-catalytic` |
| other | everything else |

**Reported per class.** For each class, the reported counts use the run's **global** top-K at 0.1 %
(no re-ranking):
- n regulators and n universe pairs;
- selected pairs and positives;
- a, OR and Fisher p.

The class assignment is label-free and lives in `data/processed/regulator_classes.tsv`, built by
`scripts/data/build_regulator_classes.py`.

### 11.2 Paralog-pair exclusion (sensitivity subset `paralog_exclusion`)

**Source.** Ensembl release 116, the release of `gene_annotation.tsv`. The intended route was BioMart
(`hsapiens_paralog_ensembl_gene`), but the query did not complete in about 40 min, so the release-pinned
Compara 116 homology dump for Homo sapiens was used instead. It is the same gene-tree homology set, and its
md5 matches Ensembl MD5SUM.
- Kept: human–human rows with `homology_type` ∈ {within_species_paralog, other_paralog}, which is
  BioMart's paralogue definition; `gene_split` is excluded.
- Count: 141,164 human paralogue pairs (23,488 within-species, 236,978 other, as directed DepMap key
  pairs).
- `other_paralog` includes distant family members (e.g. all receptor tyrosine kinases for ERBB2), so the
  exclusion is broad. The table keeps `homology_type` and % identity, so a stricter within-species-only
  variant can be computed from it.

The data are downloaded by `scripts/data/download_paralogs.py`; the checksum is in DATA.md and
`annotation_manifest.tsv`.

**Mapping.** Gene key → `ensembl_gene_id` (`gene_annotation.tsv`). An (A, B) pair is a paralog pair if
either direction is listed. `data/processed/paralog_pairs.tsv` holds the symmetric list of gene-key pairs.

**Exclusion.** Universe pairs (A, B) that are paralog pairs are removed from the universe:
- their scores are masked like same-chromosome pairs, and their positives are dropped;
- the top-K is re-ranked within the reduced universe;
- N and every structural count (per-arm × bin, per-stratum, and per-permutation pair sets) are corrected
  exactly. A brute-force test checks the permutation pair set.

### 11.3 Final-run selection and summary tables (D8)

**Selection file.** `scripts/final_eval.py` reads `results/final/selected.yaml`; its schema is in the
script docstring. It lists:
- `selected`;
- ≤ 3 `ablations` and the `reported` variants, with display `labels`;
- `baselines`;
- `controls`: score runs such as NC1 and C2, plus the built-ins PC2 and NC2 (`nc2_run`, a run without G2);
- `c1`: the runs that get the C1 pair-shuffle control;
- `sensitivity` runs;
- `subsets` (`paralog_exclusion`, `class_i|ii|iii`, `nc2`, `no_split_sibling`, `hub_excluded`) and the
  `subset_runs` they apply to;
- the `lines`.

**Subset handling.**
- `hub_excluded` is reported from every run's own metrics (`OR_hub_excluded`); it needs no extra
  evaluation.
- `no_split_sibling` always runs for the selected method.
- Duplicate run ids are evaluated once.

**Draft.** `results/final/selected.DRAFT.yaml` is the D7/D8 draft. The PI copies it to `selected.yaml` and
commits it. `--selected` is accepted only with `--mock`, so the dry run can use the draft.

**Outputs** (in `results/final/`, or `results/final_mock/` for the mock):
- `summary.json`;
- `summary.tsv` and `summary.md` (`depmaptr.report`): one row per role × run × line × subset × metric
  set.

| metric set | contents |
|---|---|
| `pooled` | ORs at the three K, Fisher p, both permutation p values and OR/median-null ratios, jackknife CI |
| `cmh` | frozen CMHs, the D6 CMH, mean AUROC, top-10 share, hub-excluded OR |
| `sign` | sign metrics |
| `callset` | call-set OR, if the score file has a call array |
| `class:<c>` | per-class table (§11.1) |
| `paired` | selected vs each run: Δlog OR, paired jackknife and bootstrap, ΔCMH |
| `C1` | pair-specific gain |
| `runtime` | evaluation wall time, plus the score run's wall time and threads from `meta.json` |
