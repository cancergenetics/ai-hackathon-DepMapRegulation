# Round 2 review: AI-agent failure modes (DESIGN v2)

## 1. Round-1 items: status

Adequately handled, no change needed: arm-level split, full-universe scoring with seeded ties, float32, sign-convention invariants, OLS known-answer test, manifests, "stop and report" rule, seeds, compute guards, label-use ban (§8.1).

Handled on paper but not enforceable as written:

**R1. [MUST] TEST labels are still readable.** §6 splits the processed labels, but the raw `data/external/pie_x_atlas_orion/de/*.parquet` (5.2 GB) still holds every TEST row, in the same working tree the sub-agents use. LOCK is a plain file anyone can delete ("removed with a logged reason" is on the honour system).
- Move the raw DE tables and `bench_test.parquet` to a sealed directory outside the project tree, readable only by the PI shell. Alternatively, encrypt them with a key held in the PI environment.
- Build `bench_train.parquet` from the sealed source by script, and delete the raw tables from agent-visible paths.
- Commit the LOCK (it lives in `results/final/`, which is not git-ignored). `make check` fails if LOCK has disappeared from git history without a decision note.

**R2. [MUST] Commit the selection before the final run.** Write `results/final/selected.yaml` (variant id, params hash, tuning-log row, git hash) and commit it before `final_eval.py` runs. `final_eval.py` refuses if the file is absent or if the variant has no tuning-log row. Otherwise "the selected method" can be chosen after peeking at any partial TEST output.

**R3. [SHOULD] Registry and resume block** (round-1 #15, #17) are marked accepted but are not in §8. Add `results/registry.tsv` (result id, script, git hash, input checksums, reviewer) and the RESEARCH_TREE "RESUME HERE" block. Also restore the blocking PI md5 check of the 26Q3 files before REPORT.md.

## 2. New failure modes from v2 additions

**N1. [MUST] The variant cap is ambiguous and the search order is undefined.**
- §5.5 lists components with parameters (G1 3×2, G2 2×2, A3 2×s×3 weights, arm handling 4+, and more). Crossing them gives far more than 30 configurations.
- Define the cap as distinct configurations evaluated, and fix the procedure: greedy forward selection in the §5.5 order (one component family at a time, best-so-far carried forward), at most N configurations per family, total ≤ 30.
- `variants.yaml` lists the family order and grids. The tuner rejects any configuration outside it.
- Define "gain exceeds the TRAIN bootstrap CI half-width" as the paired-bootstrap CI of the log-OR difference against the incumbent (a half-width of the OR's own CI is not a paired quantity). The rule is applied once per family.

**N2. [MUST] The permutation test is not yet specified well enough to be correct.**
- "Pairs that become cis are counted as negatives" removes positives and adds noise. Permute within chromosome-compatible strata, or drop incompatible pairs and keep the positive count.
- The per-regulator score histograms (§10) are built on true labels and cannot evaluate permuted labels. Specify lookup from the dense float32 matrix.
- Calibrate the test itself: run it on B0 and NC1 scores over ≥ 200 seeds and check that the false-positive rate at 0.05 is ≈ 5% (unit test on SIM).
- Report p as (1+r)/(1001), and never claim p below that.

**N3. [SHOULD] The 2-fold CV inside TRAIN adds noise, not protection.**
- No parameter is fitted to labels, so CV only splits the OR, and each fold has about one-third of regulators. A single hub knockdown (up to 5,000 hits) can dominate a fold.
- Either drop it or keep it with: fixed fold assignment hash-locked in FROZEN; geometric mean of fold ORs; a requirement that the winner beat the incumbent in both folds; and the top-10-regulator share reported for each fold.
- The CMH non-inferiority "not lower than M1" needs a CI: non-inferior if the paired-bootstrap lower bound of the log-ratio exceeds −δ, with δ pre-declared.

**N4. [MUST] G2 gain can be a benchmark artefact, not biology.**
- Own-gene non-expression drives both G2 and the benchmark: non-expressed regulators have a 0.2% hit rate against 0.45% for expressed ones, and CRISPRi cannot knock down a gene that is not expressed.
- G2 will raise the primary pooled OR by re-weighting the universe composition, and it can make NC2 look like OR ≈ 1 by construction.
- Required: add own-gene status (i/ii/iii) as a CMH stratum, and apply the non-inferiority guard within class (i). The selection rule for G2 uses the within-class (i) pooled OR.
- Hard-coded IDs (ACH-000971, ACH-001085) must be asserted against `Model.csv` names (HCT116, HEK293); HEK293 is not HEK293T, so state this proxy in the report.
- Label G2 as a context input, and do not include it in "headline method, context-free" statements.

**N5. [SHOULD] NaN → minimum score can hide bugs.** Assert a per-chromosome NaN fraction below a pinned threshold (e.g. 1%) in the manifest. A bug that NaNs a whole chromosome would otherwise silently become "no signal".

**N6. [SHOULD] Pinned manifest ranges.** Ranges for n, regulators and targets cannot be pinned before T3 runs. Freeze them in `expected_shapes.yaml` right after T3, and include that file in FROZEN.

**N7. [SHOULD] Independent reviewer needs independent metric code.** A reviewer reading the same `benchmark.py` repeats its bugs. Require a 30-line reviewer-written top-K OR and sign calculation from the raw TRAIN table, compared with the harness output on the real data. Define "clean checkout" as the git tree plus verified data symlinks, since `data/` is not committed.

**N8. [SHOULD] T1 (TP53 calls).** A new download by an agent: require DATA.md row with checksum, licence and release. State the missingness for models newer than the file (indicator is already in §4); assert known controls (HCT116 WT, HEK293 and TP53-null lines as checks).

**N9. [COULD] Tuning-log completeness.** `evaluate()` must have no non-logging path. Log crashes and repeats of the same configuration, and report both distinct configurations and total calls. `make check` compares log rows with variants.yaml.
