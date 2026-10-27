# Round 1 review: AI-agent failure modes

Scope: how agent execution can break the design, and whether the §8 guards are enforceable. Statistical content of §5 is left to other reviewers. Overall: §8 states sound rules but enforces them by convention. "The evaluation function refuses TEST unless `--final`" is a flag an agent can pass in good faith ("just a sanity check"). Each change below makes a guard mechanical and checkable in one command.

## A. Leakage and peeking

**1. [MUST] Physically separate TEST labels.** At freeze, a one-off script writes `bench_train.parquet` and `bench_test.parquet` (partitioned by regulator). `benchmark.py` exposes `load_labels("train")` only; the TEST file is read only by a PI-run `make final`, which writes `results/final/LOCK.json` (method id, git hash, split hash) and refuses a second run without a logged override. A data file cannot be passed by accident; a flag can. DATA_PROFILE was written from the full tables, so any later "quick look" at HCT116 hits would include TEST rows.

**2. [MUST] Hash-lock frozen artefacts.** Commit `results/FROZEN.sha256` (split file, `benchmark.py`, universe table, variant list). `tests/test_frozen.py` recomputes hashes and every evaluation entry point runs it. Changing a hash needs a decision-note ID, logged.

**3. [MUST] Extend §8.1 to target-level and post-hoc label use.** "Never as features" does not stop an agent dropping targets that are "noisy hubs in TRAIN" or calibrating per target on TRAIN positives. Targets are shared between TRAIN and TEST regulators (10% of regulators give >90% of hits), so hub status carries over: an indirect leak. Rule: tuning may select only global scalars or categoricals from the registered list; no per-gene, per-target or per-regulator quantity may be derived from benchmark labels. Gene lists (stress or proliferation signatures) must cite a non-Perturb-seq source in code, never agent recall of X-Atlas results.

**4. [MUST] Split by chromosome arm, not 5 Mb bins.** CN variance is dominated by arm-level events, so a TEST regulator at 8q has nearly the same CN vector, hence the same DepMap trans profile, as TRAIN regulators in neighbouring 8q bins. Tuning on those tunes on TEST's signal. Correlated arms (1q+/16q−) tie arms across chromosomes too. Assign TEST per arm, stratified by pair count, and report max |r(CN)| between TEST and TRAIN regulators as a check. Bootstrap CIs must resample these same units.

**5. [SHOULD] Pre-register variants and the selection rule.** Commit `results/variants.yaml` (id, params, hard cap ~30) plus the selection rule (best TRAIN primary OR; ties go to the simpler or cheaper variant). `tune.py` refuses unregistered ids; additions need a dated amendment. §5 lists many interacting options (S1/S2, transform, GLS, tau_cis, stratification, weighting, attribution a/b), which is an unreported grid search if uncapped.

**6. [SHOULD] Automatic tuning log.** `evaluate()` itself appends to `results/tuning_log.tsv` (variant, params hash, git hash, data checksums, seed, TRAIN metrics, wall time). A number with no log row does not count. REPORT.md derives the variant count from the log.

**7. [COULD] HCT116 caveats.** HCT116 is itself a DepMap model: not label leakage, but not independent. Run one sensitivity analysis without it, and label any HCT116-like weighting as benchmark-specific, outside the headline method.

## B. Metric gaming

**8. [MUST] Score the full universe for every method.** Gating (§5.2) or exclusions could otherwise shrink the universe or K. Unscored pairs receive the minimum score; universe and K are fixed; ties are broken by a seeded permutation; the harness asserts `n_scored == n_universe` or reports the missing count.

**9. [MUST] Store and rank by float32/64 |t|, not float16 or p.** §10 stores float16 scores: about 3 significant digits, so many ties at the K boundary (K = 0.1% of ~60M pairs is ~60k) and corrupted top-K. p-values underflow to 0 for strong effects. Harness reports the tie count at the cutoff and fails above a small threshold.

**10. [MUST] Honest uncertainty for the primary metric.** Pairs sharing a regulator or a target are not independent, so the Fisher p is hugely inflated. Headline significance should come from an arm-block permutation or bootstrap; Fisher p is descriptive. Declare in §6 that K = 0.1% of the TEST-universe count is primary and identical for all methods and controls. Print the regulator-stratified CMH OR and a "share of top-K positives from the top 3 regulators" column next to every OR.

**11. [SHOULD] Definitions as code constants with size assertions.** `P_THRESHOLD`, `LFC_MIN`, `CONTROL_MEAN_MIN` and the secondary filter are constants in `benchmark.py`. A test asserts universe size and positive count equal fixed numbers (cf. DATA_PROFILE: 62.4M pairs and 131,245 positives for the ≥20-line preview; the final definition must be recomputed and frozen).

## C. Silent data bugs

**12. [MUST] Real-data invariants asserted in code at pipeline start.**
- Sign convention end to end: the median LFC of the knocked-down gene's own feature (significant, expressed regulators) is negative; the median cis r(CN_A, y_A) is positive (PC1). A convention error flips one of them.
- IDs: 1:1 mapping after Entrez/Ensembl/HGNC rescue, no duplicated mapped IDs, logged counts of dropped and ambiguous genes. The 561 bare `ENSG` CN columns and the `_PAR_Y` ones are explicitly mapped or dropped.
- Row order of CN, expression and covariates is asserted equal by model ID.
- No NaN or inf in X, Y or Z after preprocessing; every target has SD > 0 after residualisation.
- Arm assignment on known genes (TP53 17p, MYC 8q, ERBB2 17q, a centromere-adjacent gene); cytoBand coordinates are 0-based.
- Every universe pair has chrom(A) ≠ chrom(B), using one annotation for both datasets.

**13. [MUST] Known-answer tests for the core algebra.** `assoc.py` closed form versus brute-force per-pair OLS (lstsq on Z, arms∉c, x_B) for ≥ 200 random pairs, agreeing to 1e-6 in t, including pairs where B's arm is near-collinear with x_B. Add SIM tests: planted effects recovered; a lineage/arm-confounded regulator null after adjustment but not under B1; sign flip inverts concordance; NC1 gives OR ≈ 1 and uniform p.

## D. Simplification, hallucination, over-claiming

**14. [MUST] Run manifests against expected shapes.** Each stage writes `manifest.json` (n_models, n_regulators, n_targets, n_pairs, covariate names, transform, variant id, git hash, seed, runtime). The PI and reviewer compare it with a frozen `expected_shapes.yaml`. This catches silent drops of LOCO arms, gene subsampling, S2 in place of S1, or Pearson in place of partial correlation to save memory. Add to every delegation brief: "If you cannot do X as specified, stop and report; do not substitute."

**15. [SHOULD] Reviewer checklist and result registry.** The reviewer re-runs from a clean checkout and compares output hashes; hand-verifies at least 3 raw pairs (sign, universe membership, label); checks the manifest against expected shapes and the tuning-log row count against the report. `results/registry.tsv` maps result id, script, git hash, input checksums and review status; REPORT.md may cite only reviewed rows.

**16. [SHOULD] Claims hygiene.** The 26Q3 identification rests on inference from markers and forum posts; make the PI md5 check a blocking item before REPORT.md. Tag consequential claims in LITERATURE.md as verified (URL fetched) or unverified, and have a second agent spot-check ~10 citations. REPORT wording: "association consistent with" unless the cis gate and attribution pass; the total-effect reading needs CN to be unconfounded after adjustment, which is only tested via SIM and NC1.

## E. Context, reproducibility, compute

**17. [SHOULD] Resume block.** A fixed "RESUME HERE" section atop RESEARCH_TREE.md: current milestone, frozen hashes, last registry ids, open decisions, do-not-repeat list. Commit and tag at each milestone (`Mx`).

**18. [SHOULD] Script-only rule.** Reported numbers come only from `scripts/` entry points writing to `results/` with a registry row. A `make check` runs tests plus the frozen-hash check before any evaluation. Commit small results after every run; `results/` is empty at present.

**19. [SHOULD] Compute guards.** Assert the estimated per-chunk memory is below budget (e.g. 20 of 32 GB) and take a wall-clock limit argument. GLS and all-pairs variants run first on a chromosome subsample, flagged in the manifest and excluded from reported metrics. At most one BLAS-heavy job at a time, with `OMP_NUM_THREADS` set.

**20. [COULD] Seeds.** One seed-list config file. Any stochastic step without a logged seed fails review.

**No change needed:** organoid exclusion, benchmark pinned to an HF commit, the no-benchmark-derived-features principle, the baseline/control set B0-B4/NC/PC/SIM, and the checksum table in DATA.md.
