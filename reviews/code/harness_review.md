# Independent review: benchmark harness (M3)

Reviewer: independent agent. It did not write the harness. Date: 2026-10-27.

**Scope**
- `src/depmaptr/benchmark.py`
- `scripts/benchmark/build_benchmark.py`
- `scripts/benchmark/seal_test.sh`
- `scripts/evaluate.py`
- `scripts/final_eval.py`
- `tests/test_benchmark*.py`, `tests/test_final_eval.py`, `tests/test_frozen.py`

**Spec:** DESIGN v3.1 §6–§8, D4, `docs/INTERFACES.md` and `docs/BENCHMARK.md`.

**Rules followed**
- No `*_TEST*` file was opened.
- Raw PIE DE was read only through polars lazy scans, filtered to TRAIN regulators or aggregated per feature.
- No real model score file was evaluated, and nothing was written to `results/`.

**Independent code** (session scratchpad, not in the repo):
- `rederive.py`
- `split_check.py`
- `metrics_check.py`
- `boot_cmh_check.py`

The existing suite passes: 34 passed, 2 skipped (`test_frozen`, not frozen yet).

## Verdict

**No blocker.** The tables, split and core metrics reproduce exactly from my own code. What remains is procedural: the seal is still pending, there is a crash hole in "run once", and the tuning log can be bypassed. There are also stale docs and some small gaps.

## What was verified (all exact matches)

### 1. Universe and positives from raw PIE (both lines)

The full sets were checked, not only a sample.

| | HCT116 | HEK293T |
|---|---|---|
| Regulators | 5,623, identical set and `pie_id`. Excluded: 61 unmapped, 16 Y/PAR/MT/unplaced, 3,287 not variable-CN | 2,688 |
| Targets | 10,209, identical set. Excluded: 728 unmapped, 24 Y/PAR/MT, 139 not measured | 11,000 |
| Universe pairs | 54,357,451 | 27,984,730 |
| X regulators | 360 | 180 |

Further checks:
- **Consistently tested:** `ref_mean` is constant within a feature in the raw DE. Every feature with `ref_mean ≥ 0.05` (HCT116 11,100; HEK293T 12,247) has exactly one row per perturbation.
- **Arm and chromosome:** chrom and arm match `frozen_flags`, and `frozen_flags` chrom matches `gene_annotation` for every regulator.
- **TRAIN positives:**
  - HCT116 has 81,258 in mine and 81,258 in the harness: 0 left-only, 0 right-only.
  - `p_adj` and `ueffect` have max abs diff 0. `lfc_sign`, `pos_primary`, `pos_ueffect` (|U| ≥ 0.2, D4) and `pos_p01` have 0 mismatches.
  - The cis exclusion removed 6,325 significant rows.
  - HEK293T has 40,343 = 40,343, also exact.
  - No positive has LFC = 0.
- **Per-regulator spot checks:** 28 HCT116 TRAIN regulators were checked one by one, and the positive sets were identical for each:
  - top hubs: POLR1G 2,285 hits, EXOSC6, DDX27, EXOSC8, DDX41;
  - X regulators: CCNQ, RBMX2 1,780, MECP2, PAK3, PLP1, GPR101, TMEM35A, USP26;
  - own class i, ii and iii, and zero-hit regulators.

  Every universe target was tested for each of them; for example, a chr X regulator has 9,904 = 10,209 − 305 trans targets.
- **Own-gene class:** all TRAIN regulators agree, with 0 mismatches. HCT116 has 1,388 / 171 / 2,247 in classes i / ii / iii; HEK293T has 147 / 38 / 1,705. No class-(iii) regulator has its own gene tested under another PIE id with the same HGNC id. The median own LFC is −27.6 on TRAIN.
- **Strata columns:** `cn_stat`/`cn_tercile`, TRAIN `degree_decile`, `expr_mean_s2`/`expr_decile` (recomputed from `depmap_expr_matched` × `panel_s2`), 5 Mb bins and `expr_context` all match. The only uncoordinated target is PAXX, as documented.
- **Manifest:** the `manifest.json` sha256s of every non-TEST output and input match the files on disk.

### 2. Split (DESIGN §6)

My reimplementation gives the same 42 arm units, the same order, triplets and draws, and the same TEST arms: 17q, 7q, 14q, 6q, 3p, Xp, 7p, 22q, 19p, 4q, 5p, 12q, 18p, 21q.
- The sort is descending by HCT116 count, ties in genomic order. The only tie is 11p and 22q, at 106 each.
- The draws use `default_rng(20261027).integers(3)`.
- The TRAIN folds alternate 1/2 in sorted order, and `train_folds.tsv` and the per-regulator `fold` agree.
- No gene is in both splits. Gene-level split membership is consistent across both lines and `regulator_split.tsv`.
- Every TRAIN-positive regulator is TRAIN.
- TEST rows of `universe_regulators_*.tsv` have empty `own_class` and `degree_decile`.

### 3. Metrics (synthetic score file on the real HCT116 TRAIN universe)

**Synthetic file**
- Rows and columns are shuffled.
- 7 regulators and 11 targets are missing, and there are 2 extra regulators.
- 2% of scores are NaN and 0.1% are exact zeros.
- Same-chromosome entries are set to 1e6.
- s_cis is random ±1.
- 40% of positives carry a planted signal that is 70% sign-concordant.
- 200 extra pairs are set exactly to the K-th |score|, giving 201 tied values at the boundary with 1 slot to fill.

**Independent implementation**
- Dense masks.
- My own blake2b + splitmix64 tie-break, written from the BENCHMARK.md text.
- scipy for the sample OR, statsmodels `StratifiedTable` for CMH.

**Results**

| Quantity | Result |
|---|---|
| Top-K selected set | identical, ties included |
| 2×2 at 0.1% | 17231 / 19434 / 64027 / 36564006 in both |
| OR | 506.336… in both, bit-identical |
| OR at 0.01% / 1% | equal |
| Fold log-ORs | equal |
| Sign-concordant OR | equal |
| Concordance / c₀ / n_overlap | equal |
| CMH (regulator × expression decile, and own class × expression decile), CMH p | equal |
| Missing-gene counts | 7 / 11 |

**Sign convention:** −sign(score)·s_cis is confirmed. NaN or zero scores give a predicted sign of 0, which is never concordant. A missing regulator gets s_cis = +1 and a NaN score.

**Bootstrap check:** `bootstrap_log_cmh`, which drives the §8.2 CMH guard and has no test in the suite, and `bootstrap_log_or` were compared with a brute-force weighted computation on 6 toy universes × 30 replicates. Max abs diff: 4e-16.

**Permutation formula:** I checked the per-replicate permutation algebra by hand (counts of pairs trans to both A and π(A), and the observed OR on the same pair set). It is correct, and the suite also checks it against brute force.

---

## Findings

### Major

**MA-1. The seal and freeze have not happened, yet model score files already exist.**

*Evidence*
- `results/FROZEN.sha256` does not exist and `~/DepMapTR_sealed` does not exist.
- These files are in the tree and readable:
  - `data/processed/bench/positives_*_TEST.parquet`;
  - `regulator_labels_*_TEST.tsv`;
  - `data/processed/pie/{tested,loose,targets,target_summary}_*`;
  - `de/*.parquet`.
- Meanwhile `data/results/scores/` holds 19 runs, including composite variants:
  - `g1_t30_w0_g2_thr1_w0`: its meta variant is `g1_t30_w0+g2_thr1_w0`, git `21933da…-dirty`;
  - `a3_r095_s2_w05_g2_thr1_w0`.
- There is no `results/tuning_log.tsv` and no `results/eval/`, so there is no record that any score was evaluated.

*Why it matters:* DESIGN §6 says the split sha256 is committed and TEST is sealed "before any model is evaluated". Right now, any process can `pd.read_parquet` the TEST labels; the `read_positives` name guard is advisory.

*Action:*
1. Run `seal_test.sh` and commit `FROZEN.sha256` before the first `evaluate.py` call.
2. Ask the model agent to confirm that its composite variants were chosen without benchmark labels (smoke tests only).

**MA-2. The "runs once" guarantee in `final_eval.py` has a crash hole, and the mock does not exercise the real TEST path.**

*How the hole arises*
- `LOCK` is written and committed only after every evaluation and comparison has finished.
- Each `bm.evaluate(..., "TEST")` call meanwhile writes `results/final/<run>/<line>_TEST.json` and appends a TEST row to `tuning_log.tsv`.
- Any exception after the first TEST evaluation (a bad subset entry, OOM in a `compare` of two 1.7 GB score files, a KeyError) therefore leaves TEST metrics visible with no LOCK.
- `check_guards` only looks for `LOCK`. A rerun with an edited and recommitted `selected.yaml` is then allowed.

*Why the mock does not catch this*
- `--mock` sets `split = "TRAIN"` inside `load_bench` and skips the TEST branch entirely. That branch covers:
  - the sealed-dir `read_positives(token)`;
  - the `regulator_labels` merge with `validate="1:1"`;
  - `set_readable`;
  - the hash check.
- That branch has run only on the toy in `test_final_eval.py`, so the first execution on real data is the one real run. BENCHMARK.md §5 says the mock "runs the same code path", which is not quite true.

*Fix:*
1. Commit a `results/final/STARTED` marker (or the LOCK itself) before unsealing.
2. Refuse to run if `results/final/*/` JSONs or TEST rows in `tuning_log.tsv` exist.
3. Make `--mock` write fold-2 TRAIN labels into a fake sealed dir under `*_TEST` names, with a fake FROZEN entry, so that it goes through the token/TEST branch.

### Minor

**MI-1. BENCHMARK.md is stale with respect to D4.** The code (`UEFFECT_MIN = 0.2`, `abs(ueffect) >= 0.2`), the tables and the committed `benchmark_summary.json` all use |Ueffect| ≥ 0.2:

| Split | `pos_ueffect` count |
|---|---:|
| HCT116 TRAIN | 38,342 |
| HCT116 TEST | 14,748 |
| HEK293T TRAIN | 6,232 |
| HEK293T TEST | 3,586 |

The docs disagree in three places:
- **§2 table:** the "pos. Ueffect" column still shows the literal-rule counts: 80,247 / 33,904 / 40,259 / 12,140.
- **§3:** still defines `pos_ueffect` as "|Ueffect − 0.5| ≥ 0.1".
- **§8 open issue 1:** still says "`pos_ueffect` implements the DESIGN text literally".

My recomputation gives 47.2% for HCT116 TRAIN (D4 says 47%) and 15.4% for HEK293T TRAIN.

**MI-2. The tuning log is automatic but easy to bypass.**
- Ways around it:
  - `evaluate(..., log=False)`;
  - building `Evaluation` and `compute_metrics` directly, as `calibrate_b0.py` legitimately does;
  - pointing `DEPMAPTR_ROOT` at a copy, which sends the log elsewhere.
- Any of these gives TRAIN metrics with no log row, while the §5.5 cap of 30 configurations and REPORT's "number of variants tried" rely on that log.
- `compare()` rows log an empty `variant` and `params_hash`.
- The literal §8.3 requirement is met.
- Suggestions:
  - restrict `log=False` to tests (an environment flag);
  - log in `compute_metrics`;
  - put the candidate and incumbent variant and params hash on compare rows.

**MI-3. The guards are speed bumps, not barriers.**
- `final_eval_token` checks only `sys.argv[0]`; the suite itself spoofs it with `monkeypatch`.
- Sealed files are owned by the same user, so a plain `chmod u+r` undoes the seal.
- The PIE DE is public at the pinned HF commit, so it can be downloaded again. `pie_x_atlas_orion/.cache/huggingface` is only 36 KB of metadata and holds no data.

That is acceptable, but REPORT should say plainly that the audit trail (FROZEN hashes, LOCK, tuning log) is the actual safeguard.

**MI-4. `final_eval.py` does not check in-tree drift, and the frozen set misses some inputs.**
- `final_eval.py` verifies only the sealed label hashes. It never checks the frozen in-tree hashes (`benchmark.py`, split, TRAIN tables, `variants.yaml`); only `tests/test_frozen.py` does, and only when someone runs pytest.
- `seal_test.sh` does not hash:
  - `scripts/final_eval.py` or `scripts/evaluate.py`;
  - `gene_annotation.tsv` or `pie_gene_map.tsv`, which define the universe;
  - `depmap_expr_matched.parquet`, which defines the expression deciles.
- Suggestion: have `final_eval.py` run the drift check itself and refuse on any mismatch.

**MI-5. Aggregate TEST-label information remains in the tree and in git.** This contradicts BENCHMARK.md §5's claim that "no label-derived TEST information stays in the tree". The risk is low (aggregates, plus one named hub).

| Location | TEST-derived content | Status |
|---|---|---|
| `results/split/benchmark_summary.json` | TEST positive counts, own-class counts, fraction LFC < 0, regulators with ≥ 1 positive | committed |
| `data/processed/bench/manifest.json` | own-class counts over the whole line, TEST included | |
| `docs/BENCHMARK.md` §2 | TEST rows | |
| `reviews/round2/data.md:25` | "PHB1 alone has 4,099"; PHB1 is a TEST regulator | |

**MI-6. The call-set OR cannot be computed in practice.**
- DESIGN §6 requires the "OR at the tool's call set (BH(z′) ≤ 0.05)".
- The harness computes it only if `scores.npz` has a `call` array.
- INTERFACES.md does not list `call`, and the existing score files (e.g. `M1`: regulators, targets, t, score, s_cis, cis_t) do not have one.
- The interface needs `call`, or a z′ field the harness can threshold.

**MI-7. PC2 and NC2 are not fully covered.**
- **PC2** (§7: HEK293T hits predicting HCT116 hits for shared TRAIN regulators) is not implemented anywhere: there is no score builder and no evaluation.
- **NC2** exists only as an optional `final_eval` subset. Three things about it are undocumented or missing:
  - `expr_context.fillna(0) < 1` treats the 3 CN-only regulators as not expressed;
  - on the HEK293T line, `expr_context` is HEK293, whereas §7 says "not expressed in DepMap HCT116";
  - `evaluate.py` offers no NC2 subset for TRAIN.

**MI-8. Sibling arms fall in different splits.**
- 10 chromosomes have one arm in TRAIN and the other in TEST: 3, 4, 5, 6, 12, 17, 18, 19, 21, X. Examples: Xp is TEST while Xq is TRAIN; 17q is TEST while 17p is TRAIN.
- Whole-chromosome CN events make regulators on sibling arms share CN vectors. For X, D3 notes that whole-X loss is common.
- That is the leakage mechanism DESIGN gave for rejecting 5 Mb bins, here at a weaker scale.
- This follows DESIGN; it is not a code deviation. REPORT should list it as a limitation.
- A cheap sensitivity analysis: TEST OR restricted to TEST arms whose sibling arm is also TEST (7p/7q) or absent (14q, 22q, 21q). This works at evaluation only.

**MI-9. Undocumented or loosely documented choices.** None of these is in §7's list.
1. **Fold ORs** use the full-TRAIN top-K restricted to each fold's regulators. This is documented in §4, not §7, and the DESIGN wording "log-OR gain > 0 in both folds" would also allow a per-fold re-rank.
2. **Sensitivity-set ORs** count primary positives that fail the stricter criterion as negatives, rather than excluding them.
3. **`manifest.json` inputs** omit `tested_<line>.parquet` and `vocab_*.tsv`, which `build_benchmark.py` reads to set own-gene class.
4. **`seal_test.sh` comment** says permissions are removed "on the sealed directory tree". Only the files are set to `a-rw`; directories stay 755, and the root is 700.

## Leakage summary (question 4)

- **`evaluate.py` and `evaluate()`/`compare()`:**
  - They accept TRAIN only.
  - `load_bench("TRAIN")` reads `universe_regulators_*` (TEST rows contain only label-free columns), `universe_targets_*` and `positives_*_TRAIN`.
  - No path reads a TEST label or TEST own-class/degree.
  - TEST regulator identities, arms and `cn_tercile`/`expr_context` are visible, but they are label-free.
- **Seal coverage:**
  - It covers the TEST label files, raw `de/*.parquet` and every `data/processed/pie/*` cache with per-perturbation rows.
  - It leaves `features_*.tsv` and `vocab_*` in place; both are label-free.
  - `preprocessed/meta.json` holds only the gene list.
  - No PIE extracts were found in other DepMapTR session scratchpads.
  - Not covered: the aggregates in MI-5, and the public re-download route in MI-3.
- **Tuning log:** it is automatic by default inside `evaluate()` and `compare()`, but it can be bypassed (MI-2).
