# DESIGN — DepMap CN → trans-expression effects (DepMapTR)

Version: **v3.1, FROZEN** (2026-10-27; after round 3: one fresh reviewer, reviews/round3/). v1 was reviewed by 7 agents (reviews/round1/) and v2 by 5 (reviews/round2/).
Every proposal and its disposition is in §13. Inputs: PROMPT.md (fixed goal), LITERATURE.md, DATA.md,
docs/DATA_PROFILE.md.

## 1. Goal and estimand

For a regulator A and a target B on different chromosomes, the tool estimates the effect of A's copy
number on B's mRNA, β_AB, together with a signed score, a test statistic and a ranking. The target is the
**total** (direct + indirect) consequence of a change in A's dosage. Estimation holds fixed B's cis CN,
the arm-level CN of other chromosomes and cell-intrinsic context.

Two caveats are stated explicitly:
- In the default (LOCO) model, β_AB is the effect of "A's dosage including the arm-level dosage it shares
  with its neighbours". Every gene on A's arm therefore inherits a similar trans signature. Attributing
  that signature to the gene actually responsible is a separate step (§5.4), and a focal-only variant
  exists (§5.5).
- The CRISPRi benchmark lowers A. The predicted knockdown direction is therefore
  **sign_kd(A→B) = −sign(β_AB) · s_A**. Here s_A = −1 only if A's own cis CN→mRNA t ≤ −3, and +1
  otherwise. About 95% of genes scale positively, so a weak cis estimate must not flip the sign at random
  (round-2 domain #1, statistics #6). t_cis is the G1 cis t (§5.3). If A has no expression column,
  t_cis = 0. For the hinge variant, β_AB is the loss slope.

## 2. Data (see DATA.md, docs/DATA_PROFILE.md)

- DepMap Public 26Q3 (inferred release, checksums in DATA.md): `OmicsCNGeneWGS` (linear relative CN,
  ploidy-normalised) and `OmicsExpressionTPMLogp1HumanProteinCodingGenes` (log2(TPM+1)). One
  `IsDefaultEntryForModel` profile per model.
- Covariates:
  - Model.csv: lineage, primary disease, sex, GrowthPattern, EngineeredModel, PatientID.
  - ModelCondition: media.
  - OmicsProfiles: RNA `Stranded`, platform.
  - 24Q4 PureCN: ploidy, WGD (60% coverage). 24Q4 MSIScore (76%).
  - TP53 functional status from the DepMap **24Q4** figshare mutation matrices (pinned; checksums in
    DATA.md). This is data task T1.
  - Arm-level CN derived from gene CN.
- Benchmark: PIE X-Atlas/Orion `de/*.parquet` (HF commit abfa3da), HCT116 and HEK293T. These are
  pre-computed perturbation–mRNA associations only.

## 3. Samples and preprocessing

- **Primary sample set S2:** all DepMap *cell lines* (no organoids) with matched default CN and RNA,
  with two exclusions:
  - EngineeredModel lines (drug-adapted, knockouts, clones);
  - all but one model per PatientID (prefer a model with PureCN, then the default).

  This gives n = 1,436 (round-2 data check); the exact n is logged in the manifest.
  Rationale: S1 (PureCN-complete, n ≈ 1,037) loses 31% of samples, drops HCT116 itself, and drops 45% of
  MSI-H lines (data review E2). **S1 is a sensitivity analysis only.**
- **CN transform:** x = log2(clip(relCN, 2^-2, 2^3)). The floor caps the leverage of homozygous deletions
  at −2. Homozygous deletion is handled explicitly through the cis basis below. Pre-listed variants:
  floor 2^-5, and rank-based inverse-normal (RINT).
- **Expression:** log2(TPM+1) as provided. Variant: RINT.
- **Frozen universe flags** (computed once on the fixed panel P = S2, before tuning, committed and
  hash-locked; no tunable may change them):
  - regulator **variable CN**: ≥ 20 lines in P with |log2 relCN| ≥ 1. For X genes the test is applied to
    **sex-centred** CN (log2 relCN minus that gene's median within the line's sex), because male X sits at
    log2 ≈ −1 (round-2 statistics #2, data #2, external PI);
  - DepMap target **measured**: present in both the CN and expression matrices, finite, SD > 0 in P.
    Targets without CN are dropped, because the brief requires each gene's cis CN.
- Chromosomes:
  - Autosomes and X are analysed. X regulators are kept, with sex as a covariate (external-PI review #12).
  - Y and PAR genes are excluded as regulators and as targets, which avoids the X/Y PAR column ambiguity
    (data review E4).
  - MT genes are excluded.

## 4. Covariates and justification

Shared covariates Z, present in every model:

| covariate | encoding | why |
|---|---|---|
| lineage | OncotreeLineage; heterogeneous lineages (Lymphoid, Myeloid, CNS/Brain, Esophagus/Stomach, Soft Tissue, Bone) split by OncotreePrimaryDisease when the level has ≥ 15 lines; levels < 15 merged into the parent lineage or "other" | the strongest expression confounder; lineage-specific arm events (CRC 18q−/20q+) create cross-chromosome CN–expression correlation |
| ploidy, WGD | PureCN values; missing → lineage-median imputation + missingness indicator | ploidy-normalised CN hides absolute dosage; WGD changes global transcription and cell size |
| genome burden **LOCO** | fraction of **autosomal** arms altered (|arm log2| > 0.3) and mean |arm log2|, computed **excluding A's chromosome**, inside the per-chromosome step | the generic aneuploidy-stress axis. The PureCN CIN/LoH values include A's own chromosome and would partly adjust away x_A (domain review #3), so they are not used |
| MSI | MSIScore > 20 flag + missingness indicator | MSI lines are near-diploid with a distinct programme; HCT116 is MSI-H |
| TP53 status | DepMap 24Q4 figshare mutation matrices: TP53 altered iff Damaging > 0 or Hotspot > 0; missingness indicator (covers 1,255/1,436 lines) | drives CN burden and p53-target expression, a major generic responder set in TP53-WT HCT116. A mutation call is upstream of expression, unlike an expression-derived p53 score (not used) |
| sex | from Model.csv; missing → indicator | X-target dosage; X regulators |
| media | base medium (RPMI/DMEM/F12/…, levels ≥ 15) + serum-free flag | culture-driven programmes |
| growth pattern | adherent/suspension/mixed/spheroid | suspension vs adherent programmes |
| RNA technical | `Stranded` flag, sequencing platform (levels ≥ 15) | known library/batch axis; confounded with sample set |
| purity | not available; cell lines assumed ≈ 1 (documented assumption) | — |

Pair-specific covariates:
- **Cis basis for B:** W_B = [x_B, 1(relCN_B < 0.25)]. This covers both the linear cis dosage and
  homozygous deletion of B, as the brief requires (statistics review #2). The deletion indicator is
  included only when ≥ 3 lines carry it. k_B and df are therefore per target (round-2 statistics #1).
- **Arm CN definition:** arm_k = median of x over the genes on arm k (UCSC GRCh38 band at the gene
  midpoint; arms with ≥ 20 genes), computed with the active CN transform. The burden thresholds apply to
  this quantity.
- **Regulators without an expression column** (CN-only genes): they fail G1 (t_cis = 0, s_A = +1) and
  count as "not expressed" for G2, and B2 is set to its minimum for them. Their number is logged in the
  manifest.
- **LOCO arm CN:** the arm-level CN of every arm *not* on A's chromosome. These absorb correlated arm events
  across chromosomes; B's own arm is included. Xp and Xq arm CN are sex-centred.

Deliberately **not** used:
- **Expression PCs or PEER factors.** They absorb genuine CN-driven trans hubs. An ablation with
  regulator-conditional factors (factors computed from Y after removing the CN effects of chromosome c)
  is pre-listed.
- **A proliferation score from expression.** Used only as a sensitivity analysis.

We also note one deliberate choice: by adjusting genome burden, the model blocks real
A → aneuploidy → B paths (e.g. CCNE1). Effects of that kind are under-estimated by design.

## 5. Statistical model

### 5.1 Core model M1 (reduced form, LOCO)
For each pair: y_B = β_AB·x_A + W_B·γ + Σ_{k∉chr(A)} δ_k·arm_k + Z·η + ε.

The fit uses Frisch–Waugh–Lovell. For each regulator chromosome c (23 groups):
- D_c = [1, Z, burden_LOCO(c), arms∉c]. Residualise Y and X on D_c using float64 pivoted QR. Log cond(D_c)
  and rank(D_c).
- Orthonormalise each target's residualised W̃_B into Q_B. The partial correlation of x_A and y_B given
  (D_c, W_B) is then
  - ⟨x̃_A, ỹ_B⟩ − (Q_B'x̃_A)'(Q_B'ỹ_B), divided by the square root of
  - (‖x̃_A‖² − ‖Q_B'x̃_A‖²)·(‖ỹ_B‖² − ‖Q_B'ỹ_B‖²).
- The cost is k+1 GEMMs per chromosome (k = 2), with df_c = n − rank(D_c) − k − 1.
- **Guards:** a pair is set to NaN where either residual norm falls below 1e-3 of its original norm, or the
  cell line count is degenerate. NaN pairs receive the minimum score in evaluation.
- t is computed from r; p-values come from log-space t tails.
- **Known-answer test:** agreement with brute-force OLS (statsmodels) to 1e-8 on ≥ 200 random real pairs.

### 5.2 Calibration and calls
- **NC1 permutation null:** permute rows of the residualised regulator matrix X̃ within lineage strata,
  jointly across genes. Keep Ỹ and W̃_B fixed (as in tensorQTL residual permutation).
- **Diagnostics:** report λ_GC and per-chromosome QQ plots for the real and NC1 scores.
- **Empirical null:** Efron central matching per regulator chromosome gives z′. BH is applied to z′ over
  the family "all analysed regulators × measured trans targets"; BY is reported as a bound. The tool's
  **call set** is BH(z′) ≤ 0.05. The number of distinct CN segments called is also reported.

### 5.3 Gates (label-free and DepMap-internal)
- **Instrument gate (G1).** The cis model is y_A ~ x_A | D_c (with A's own arm handled as in the variant).
  Regulators need a positive cis slope with first-stage F = t²_cis ≥ τ (τ ∈ {10, 30, 100}; pre-listed).
  Non-gated regulators keep their score multiplied by a weight w_G1 ∈ {0, 0.5} (pre-listed).
  - β_IV = β_AB/β_cis is reported with an Anderson–Rubin/Fieller interval, labelled "per unit A mRNA,
    if A is the sole mediator".
  - The reduced-form t is the AR test of β_IV = 0.
- **Context-expression gate (G2).** Is A expressed in the benchmark's context, according to DepMap?
  - Expression of A in DepMap HCT116 (ACH-000971) for the HCT116 benchmark, and in HEK293 (ACH-001085)
    for HEK293T.
  - Thresholds: log2(TPM+1) ∈ {1, 2}, pre-listed. Pairs from non-expressed regulators get weight
    w_G2 ∈ {0, 0.5}.
  - This uses no benchmark labels. It is a tool input ("context cell line"). G2 is reported as its own
    variant, because it makes predictions context-specific.
  - The context model IDs are asserted against Model.csv (CellLineName HCT 116 / HEK 293).
  - Own-gene benchmark status correlates with G2. Its gain is therefore also reported *within* regulator
    class (i) and in a CMH analysis that includes own-gene status as a stratum (round-2 AI #N4). A true
    non-expressed regulator cannot be knocked down, so G2 is legitimate DepMap-derived prior knowledge.
    It is still kept separately visible.

### 5.4 Co-segment attribution (A3)
1. Greedy sequential segmentation along each chromosome. A gene joins the current segment if
   r(CN, segment mean) ≥ ρ, with ρ ∈ {0.95, 0.98} pre-listed.
2. Within a segment of m genes (**all** DepMap genes in the segment, not only perturbed ones, so the
   benchmark design cannot leak in; round-2 external PI #3), compute a single-effect (SuSiE L = 1)
   posterior for each target B:
   α_{g|B} ∝ π_g · ABF(z_gB) · ABF(z_co,gB)^w, where:
   - ABF(z) = √(1−λ)·exp(λz²/2), with λ = σ₀²/(σ₀² + 1/df) and **σ₀ = 0.05 fixed a priori**;
   - the co-expression exponent w_co ∈ {0, 0.25, 0.5} (written w above);
   - if Σπ = 0 within a segment, π is uniform;
   - π_g ∝ exp(t_cis,g / s) · w_G2(g), with s ∈ {2, 5} and w_G2 the context weight of §5.3.

   Singleton segments get α = 1. L = 1 does not remove arm-wide signal shared across segments; the focal
   variant does that.
3. **Residual co-expression** z_co,gB = atanh(r)·√(df − 1), where r = corr(ỹ_B, residual of y_g on
   x_g and D_c | W_B). It
   is co-expression, not mediation, and is used **only** inside the formula above for within-segment
   redistribution, never as an additive main effect. It never conditions on y_A as a covariate for B,
   which would make y_A a collider.
4. The signed pair score is t_AB · α_{A|B}.

### 5.6 Final score
S_AB = t_AB · w_G1(A) · w_G2(A) · α_{A|B}:
- inactive components equal 1;
- when A3 is active, w_G2 enters only through π_g, not again as a multiplier;
- direction-agnostic ranking uses |S_AB|;
- the predicted knockdown sign is −sign(t_AB)·s_A (§1).

Under target standardisation, t_AB is replaced by its z-score across regulators before S is formed.

### 5.5 Pre-registered variants (`results/variants.yaml`, hash-locked)
**Cap: at most 30 configurations evaluated on TRAIN in total.** Tuning is **greedy, one tier at a time, in
the order below**. Each tier starts from the current best configuration, tries that tier's options, and
adopts one only under the rule in §8.2 (round-2 external PI #2, AI #N1).
1. Arm handling:
   - LOCO (default);
   - focal = residualise on all arms including A's own;
   - leave-one-arm-out;
   - LOCO CN-PCs (k = 20) instead of arms.
2. G1 instrument gate: τ ∈ {10, 30, 100} × w ∈ {0, 0.5}.
3. G2 context gate: threshold ∈ {1, 2} × w ∈ {0, 0.5}.
4. A3 attribution: ρ ∈ {0.95, 0.98}; s ∈ {2, 5}; co-expression weight w ∈ {0, 0.25, 0.5}. Searched as a
   one-at-a-time sweep from (0.95, 2, 0), not a full grid.
5. Target standardisation: z-score each target's t across regulators before pooled ranking.
6. Dosage shape: hinge at copy-neutral (log2 relCN = 0), with separate loss and gain slopes; the loss slope
   is used for ranking. This is motivated by CRISPRi mimicking loss.
7. Transforms: CN floor 2^-5; RINT of X and Y.

Sensitivity analyses (not tuned; reported on TRAIN and in the final run):
- sample set S1;
- leave-HCT116-out;
- regulator-conditional hidden factors;
- paralog-pair exclusion;
- CN-outlier exclusion.

GLS with a LOCO CN-kinship matrix is deferred (only if time allows).

Dropped:
- the lineage-stratified meta-analysis (underdetermined, n = 20–180 per lineage vs p ≈ 45);
- context weighting toward colorectal/MSI-H lines (those lines are CN-quiet).

## 6. Benchmark definitions (frozen before any model is run; `src/depmaptr/benchmark.py`, hash-locked)

- **Terminology.** PIE `target` = perturbed gene = our regulator A. PIE `feature` = measured gene = our
  target B. "BH within target" means BH within each perturbation.
- **Operational exclusions** from the literal universe: Y, PAR, MT and unmapped genes. Their counts are
  reported.
- **Cell lines.** HCT116 is primary: a DepMap model, colorectal MSI-H, 2× the targets, and a plausible null
  level among non-expressed targets. HEK293T is secondary: not a DepMap model, and 83% of its trans hits
  come from non-expressed targets.
- **Positive (direction-agnostic):** a trans pair with PIE `p_adj < 0.05` (BH within target, as released).
  - Sensitivity sets: |Ueffect − 0.5| ≥ 0.1 & p_adj < 0.05; p_adj < 0.01.
  - LFC magnitude is never used, because of the pseudocount artefact.
- **Sign:**
  - The benchmark sign is sign(LFC).
  - The predicted knockdown sign is −sign(β_AB)·s_A, with s_A as defined in §1.
  - A pair is **concordant** iff predicted knockdown sign = sign(LFC). This is unit-tested, including the
    case of a regulator with a negative cis slope.
- **Universe (primary, literal to the brief):** trans pairs (both genes placed, chrom(A) ≠ chrom(B))
  where:
  - A is a perturbed gene (PIE `target`) in that line's DE table and maps 1:1 to DepMap;
  - A is on an autosome or X and has a frozen variable-CN flag;
  - B is consistently tested in that line (control mean ≥ 0.05) and has a frozen DepMap-measured flag;
  - A ≠ B.

  Every method scores **every** universe pair. Missing or NaN scores get the minimum, and ties are
  broken by a seeded hash.
- **Secondary universe.** Regulators are labelled by the benchmark's own-gene status: (i) own feature
  tested and knocked down (p_adj < 0.05, LFC < 0); (ii) tested but not knocked down; (iii) not tested
  (not expressed). Results are reported for (i) and for (iii). Set (iii) is an empirical-null set where
  the expected OR is ≈ 1.
- **Held-out split** (frozen before any tuning; seed 20261027):
  - The split unit is the **chromosome arm**. Arm-level events dominate CN covariance, and 5 Mb bins would
    leak through near-identical CN vectors.
  - The arms of chr1–22 and X that contain ≥ 1 universe regulator (either cell line) are sorted by the
    number of HCT116 universe regulators. In consecutive triplets, one arm is randomly assigned to TEST,
    giving ≈ 1/3 TEST.
  - `results/split/regulator_split.tsv` (gene → arm → split) is committed with its sha256 in
    `results/FROZEN.sha256` before any model is evaluated.
  - At freeze, the benchmark labels are physically split into `data/processed/bench_train.parquet` and
    `bench_test.parquet`.
  - The test file is then **sealed**, together with the raw `de/*.parquet` and any cache holding TEST
    rows. Sealing means:
    - moving them out of the project tree to `~/DepMapTR_sealed/`;
    - removing their read permissions;
    - recording their sha256 in `results/FROZEN.sha256`.
  - Tuning code reads only `bench_train.parquet` (round-2 AI #R1).
  - `scripts/final_eval.py` refuses to run unless `results/final/selected.yaml` (the selected method and
    the ≤ 3 ablations) is committed. It restores read access, runs once, and writes and commits
    `results/final/LOCK`. Any rerun refuses while the LOCK is committed (round-2 AI #R2).
  - A dry run of `final_eval.py` on a **mock** test file (TRAIN arms relabelled) is done before the real
    run.
- **Metrics.** Each is computed per cell line and per split.
  - **Primary:** pooled Fisher OR among the top K = 0.1% of universe pairs of the evaluated split and
    cell line, ranked by |S| (HCT116 TEST).
  - **Resampling conventions** (round-3 #4):
    - In every jackknife, bootstrap and permutation replicate, the top-K score threshold is fixed at its
      full-split value.
    - Replicate ORs are computed from 2×2 counts per (regulator arm × target 5 Mb bin).
    - Permutations reassign labels only among regulators of the same split and cell line.
  - **Primary significance:** a regulator-label permutation p (1,000 permutations).
    - Benchmark label rows are reassigned among universe regulators within strata (own-gene status ×
      CN-variability tercile).
    - Each permuted OR is computed on the pairs that are trans for both A and π(A). The observed OR in
      each permutation uses that same pair set, so pairs are never recoded as negatives, and
      p = (1 + #{OR_null ≥ OR_obs})/1001.
    - This preserves the score structure and target hubness. It does **not** preserve regulator benchmark
      degree, so it credits both regulator selection and pair-specific knowledge (wording corrected,
      round-2 statistics #3).
    - A secondary, evaluation-only permutation within strata of benchmark-degree decile × own-gene status
      isolates pair-specific signal.
    - The permutation tests are calibrated on B0 and NC1 (≥ 200 seeds: the p-values must be ≈ uniform).
    - The ratio OR / median null OR is reported.
  - **CIs:** a delete-one-arm jackknife on log-OR with t_{G−1} critical values (G = number of TEST arms,
    ≈ 14), with targets held fixed. Methods are compared with the same paired jackknife. The two-way bootstrap (regulator arms × target bins) is a sensitivity analysis (round-2
    statistics #4). REPORT.md must state that TEST precision is limited by G.
  - **Fisher p:** reported, labelled "assumes independence".
  - **Secondary:**
    - OR at K ∈ {0.01%, 1%} and at the tool's call set (BH(z′) ≤ 0.05);
    - Cochran–Mantel–Haenszel OR stratified by regulator × target-expression decile, which reports
      pair-specific knowledge net of regulator selection and target detectability. A second CMH analysis
      uses own-gene status × target-expression decile strata;
    - per-regulator AUROC (mean over regulators with ≥ 10 positives);
    - share of top-K pairs from the top 10 regulators;
    - hub-excluded sensitivity: drop the top 1% of regulators and of targets by benchmark degree. This is
      evaluation-only and never used for tuning.
  - **Sign:**
    - the sign-concordant OR (positives = benchmark-positive and concordant);
    - the concordance rate among overlapping pairs vs the independence null c₀ = qu + (1−q)(1−u), where
      q = the fraction of overlapping pairs whose predicted knockdown sign is down and u = the fraction
      with LFC < 0 (unit-tested);
    - a 2×2 sign OR;
    - a p-value from within-regulator sign permutation;
    - regulator-weighted concordance (regulators with ≥ 5 overlaps).
  - **Runtime:** wall clock and threads for every method.

## 7. Baselines and controls

| id | description | purpose |
|---|---|---|
| B0 | random scores (seed list) | null, OR ≈ 1 |
| B1 | naive Pearson r(x_A, y_B), no covariates | value of the covariates |
| B2 | co-expression r(y_A, y_B), no CN | does CN add anything over co-expression? |
| B3 | regulator hub: rank by cis t_A × SD(x_A), target-agnostic | degree baseline |
| B4 | target hub: SD of y_B across DepMap | degree baseline |
| B4e | target detectability: mean y_B in DepMap HCT116 | expression-power baseline |
| B5 | rank(B3) × rank(B4e) | combined degree baseline; the method must beat it on pooled **and** CMH OR |
| B6 | arm-only: r(arm CN of A's arm, y_B) adjusted as M1 | is gene-level CN beyond the arm informative? |
| NC1 | residual-permutation null (§5.2) | negative control, OR ≈ 1 |
| NC2 | class (iii) regulators (not expressed in the benchmark) that are also not expressed in DepMap HCT116 (log2 TPM+1 < 1); evaluated **without** G2 | biological null |
| PC1 | cis pairs (own CN → own mRNA) | positive control for the CN pipeline (outside the metrics) |
| PC2 | HEK293T hits predicting HCT116 hits for shared TRAIN regulators | empirical reference OR |
| SIM | toy simulation: planted trans effects, arm co-selection across chromosomes, lineage confounders, co-segment neighbours, a negative-cis regulator, and homozygous deletions | recovery, calibration, sign convention, confounder removal |

## 8. Rules (overfitting and AI-agent failure modes)

1. Benchmark labels are used **only** to choose among the pre-listed variants (§5.5) on TRAIN regulators.
   They are never used as features, never for target- or regulator-level decisions (no dropping of
   "noisy" targets, no target-hub priors), and never for calibration. DepMap CRISPR essentiality, LINCS,
   other Perturb-seq data and perturbation-trained embeddings are excluded from the primary method.
2. **Tuning objective:** the pooled top-K OR on HCT116 TRAIN.
   - **Adoption rule** (round-2 external PI #1): a variant replaces the current best only if
     - its log-OR gain is > 0 in **both** folds of a fixed 2-fold arm-block split of TRAIN, **and**
     - the 90% lower bound of a paired two-way bootstrap (TRAIN regulator arms × target bins, 1,000
       replicates) of the log-OR difference is > 0.
   - **CMH guard:** the candidate must not be significantly worse than the current incumbent on
     regulator-stratified CMH OR (90% upper bound of the paired difference > 0).
   - The 2-fold split alternates TRAIN arms in the §6 sorted order (seed 20261027) and is committed in
     FROZEN.sha256.
   - If several options in a tier pass, the one with the largest TRAIN Δlog-OR is adopted.
3. `evaluate()` appends automatically to `results/tuning_log.tsv`: variant id, params hash, git hash,
   data checksums, seed, TRAIN metrics, wall time. REPORT.md reports the number of variants tried.
4. **One final TEST run** (`scripts/final_eval.py`). It evaluates:
   - the selected method;
   - all §7 baselines and controls;
   - three ablations relative to the selection: selected minus G2; selected minus G1; selected with the
     other arm handling (LOCO↔focal);
   - the §5.5 sensitivity analyses, listed in `selected.yaml`.

   All results are reported.
5. `results/FROZEN.sha256` covers the split, `benchmark.py`, the universe-flag tables and `variants.yaml`.
   `tests/test_frozen.py` fails on any drift. Changes after the freeze need a decision note in
   RESEARCH_TREE.md and must be reported.
6. Real-data invariants are asserted at pipeline start:
   - the median on-target LFC is negative;
   - the median cis r is positive;
   - gene IDs are 1:1 and the row orders of CN, RNA and covariates are aligned;
   - there are no NaNs after preprocessing;
   - known genes are on the right arms (ERBB2 → 17q, MYC → 8q, CDKN2A → 9p);
   - chrom(A) ≠ chrom(B) for every evaluated pair.
7. Each stage writes `manifest.json` (n_models, n_regulators, n_targets, n_pairs, covariates, variant id).
   Its expected ranges are pinned in tests. Delegation briefs say: "if something does not fit, stop and
   report; do not substitute or simplify silently."
8. Every result relied on is reviewed by an agent that did not produce it. The reviewer reruns it from a
   clean checkout, compares hashes, and hand-checks ≥ 3 raw pairs.
9. Seeds are in `config/seeds.yaml`: [20261027, 1, 2, 3, 4]. The split uses 20261027 only. Permutation
   and NC1 seeds are derived with `numpy.random.SeedSequence(20261027).spawn(n)`, and the spawn index is
   logged.
11. At each milestone: a decision note and a two-paragraph status in RESEARCH_TREE.md; prune LESSONS.md;
    ask a fresh agent "what should we stop doing / what worked well".
10. Compute guards: an estimated memory per chunk < 20 GB, assertions on that, and a wall-clock limit
    argument.

## 9. Software layout

- `src/depmaptr/`:
  - `io.py`: matched matrices + covariates, invariants;
  - `covariates.py`;
  - `assoc.py`: generic engine for pair-specific regressors, chunked, float64;
  - `calibrate.py`;
  - `gates.py`;
  - `attribution.py`;
  - `benchmark.py`: universe, labels, sparse metrics, permutation, bootstrap, split guard;
  - `simulate.py`.
- `scripts/`: CLI entry points (`run_model.py --variant`, `evaluate.py`, `final_eval.py`) and
  `make check` (tests + frozen hashes).
- `results/`: committed small outputs.
- `data/results/`: dense float32 score matrices (ignored).

## 10. Compute plan

Measured by the algorithms review:
- M1 over ≈ 2×10^8 pairs: about 1 minute (BLAS);
- 100 permutation passes: about 35 minutes;
- ≈ 30 variants plus final: 1–3 hours.

Scores are stored as dense float32 per regulator chromosome, about 1 GB in total. Evaluation uses sparse
positives plus per-regulator score histograms, by label class, for the fast OR, CMH and bootstrap
computations. numpy here uses Apple Accelerate: matmul warnings are spurious and must never become
errors.

## 11. Data tasks before modelling
- T1: TP53 status (DepMap mutation calls); latest programmatically available release, in DATA.md.
- T2: fix DATA_PROFILE §9 (687 same-WGS-profile, not 905; data review E1).
- T3: recompute the variable-CN and measured flags on panel P; build the S2 sample table
  (dedup, engineered exclusions).
- T4: RNA `Stranded`/platform covariates.
- T5: LITERATURE fixes (ERBB2 17q; CausalBench conclusion updated; unverified numbers marked).

## 12. Milestones

M3: harness + frozen split + benchmark tests + SIM → M4: baselines + M1 + NC1 → M5: variants on TRAIN →
M6: one final TEST run → M7: REPORT.md.

**Minimum viable pipeline** (built first and run end-to-end before any extras; round-2 external PI #4):
- M1 + G1 + G2;
- baselines B0–B6 and NC1;
- SIM;
- the permutation test and jackknife;
- a mock `final_eval.py` dry run.

**Deferrable** (only after the minimum viable pipeline works):
- Efron empirical null. Until it is built, the call set is {|t| ≥ t*}, where t* is the smallest threshold
  with (#NC1 pairs ≥ t* / n_perm) / (#real pairs ≥ t*) ≤ 0.05;
- AR/Fieller intervals;
- GLS;
- the S1, hidden-factor and paralog sensitivity analyses;
- the regulator-class stratification.

## 13. Review log

### Round 1 (7 reviewers; reviews/round1/*.md). Dispositions:
| proposal (reviewer) | decision |
|---|---|
| Separate train/test label files, LOCK, hash-lock (AI #1, #2) | **accepted** §6, §8.5 |
| Ban label-derived target/regulator decisions (AI #3) | **accepted** §8.1 |
| Split by arm, not 5 Mb (AI #4, ext-PI #6, prior #4, data #4) | **accepted**, arm triplet assignment §6 |
| Score full universe; float32; seeded ties (AI #8, #9; alg #4) | **accepted** §6, §10 |
| Dependence-aware significance (ext-PI #2, stats #9, prior #2, AI #10) | **accepted**: regulator-label permutation (stats #9a) as primary, rather than within-regulator score shuffling (ext-PI #2), because it preserves target hubness; two-way bootstrap |
| Sign-concordance null c₀ + within-regulator permutation (ext-PI #3, stats #10) | **accepted** |
| Knockdown sign from the IV sign (domain #7, stats #6) | **accepted** §1, §6 |
| AR/Fieller for IV, F ≥ 10 gate (stats #6) | **accepted** §5.3 |
| Context-expression gate (ext-PI #4, domain #10) | **accepted** G2 |
| Focal vs arm variants, B6 (ext-PI #5, stats #3, alg #7) | **accepted** §5.5, B6 |
| Selection rule, single TEST run, ≤ 3 ablations (ext-PI #1) | **accepted** §8.2, §8.4 |
| Tuning objective = CMH OR (stats #11) | **modified**: the brief's primary is pooled enrichment, so the pooled OR is the objective and CMH is a non-inferiority guard plus the B5 requirement |
| Gated scoring as primary (prior #1) | **modified**: G1 is first in the variant ranking; the primary is chosen on TRAIN by rule §8.2 |
| Target-detectability stratification and baseline (prior #3) | **accepted** CMH strata, B4e |
| Combined hub baseline (ext-PI #11) | **accepted** B5 |
| df/rank/missingness exact; OLS known-answer test (stats #1, AI #13) | **accepted** §5.1 |
| Cis basis with homozygous-deletion indicator (stats #2) | **accepted** §4 |
| Empirical null, λ_GC, residual-permutation NC1 (stats #5, alg #11) | **accepted** §5.2 |
| "Mediation" renamed residual co-expression; within-segment use only (stats #7, ext-PI #9, prior #6) | **accepted** §5.4 |
| SuSiE L = 1 attribution, greedy segments (stats #8, alg #8, #9) | **accepted** §5.4 |
| GLS LOCO kinship, single δ (stats #12, alg #10) | **accepted**, lowest priority |
| Lineage-stratified meta-analysis (stats #14, alg #13) | **dropped** |
| Hinge / loss-weighted slope (domain #8, ext-PI #14) | **accepted** variant |
| TP53, RNA technical covariates, LOCO burden (domain #1–3, data #5, ext-PI #10) | **accepted** §4 |
| Finer lineage (domain #5) | **accepted** §4 |
| Heterogeneity-aware pooling instead of colorectal weighting (domain #11) | **partly accepted**: colorectal weighting dropped; heterogeneity pooling not added (no budget, would need lineage-level fits) |
| S2 primary, frozen universe flags (data #1, #2, E2, E3) | **accepted** §3 |
| Dedup patients, drop engineered (data #3, ext-PI #10) | **accepted** §3 |
| CN floor leverage (data #6, stats #4, E5) | **accepted**: floor 2^-2 + cis indicator; 2^-5 variant |
| Keep X regulators (ext-PI #12); PAR handling (data E4) | **accepted** §3 |
| Ueffect sensitivity positives; three-state own-gene status; hub-excluded sensitivity (domain #14–16) | **accepted** §6 |
| Pre-registered regulator-class stratification (domain #9) | **accepted** as a REPORT analysis (HGNC gene groups: TF, RBP, kinase, ribosomal, complex member); not used for tuning |
| Goh reverse check on co-expression (domain #12) | **deferred**: co-expression already restricted to within-segment use |
| DepMap-only regulator prior (prior #5) | **covered** by B3/G1 (cis × CN SD); no separate prior |
| Paralog-pair flag (ext-PI #15) | **accepted** as a sensitivity analysis |
| Leave-HCT116-out sensitivity (AI #7, prior #9) | **accepted** as a sensitivity analysis |
| Manifests, invariants, reviewer checklist, resume block, compute guards (AI #12, #14, #15, #17, #19) | **accepted** §8 |
| Per-regulator score histograms (alg #5) | **accepted** §10 |
| Heteroskedasticity-robust t (stats #13) | **accepted** as a check on the final method only |
| TPM compositional artefact (domain #13) | **noted** in limitations; checked by excluding the top 1% of CN outliers in a sensitivity analysis |
| LITERATURE fixes (prior work) | **accepted** T5 (done) |

### Round 2 (5 reviewers; reviews/round2/*.md). Dispositions:
| proposal (reviewer) | decision |
|---|---|
| Sign rule with shrunk cis sign (domain #1, stats #6) | **accepted**: s_A = −1 only if t_cis ≤ −3 (§1) |
| LOCO burden on autosomes only (domain #2) | **accepted** §4 |
| X variable-CN on sex-centred CN; sex-centre Xp/Xq in LOCO (stats #2, data #2, ext-PI) | **accepted** §3, §4 |
| Conditional homozygous-deletion indicator; per-target k_B and df (stats #1) | **accepted** §4 |
| Permutation on pairs trans for both A and π(A); corrected interpretation; degree-stratified secondary; calibrate on B0/NC1 (stats #3, AI #N2) | **accepted** §6 |
| Jackknife over TEST arms with t_{G−1} (stats #4) | **accepted** §6 |
| SuSiE L = 1 formula with fixed σ₀ (stats #5) | **accepted** §5.4 |
| Attribution normalised over all segment genes (ext-PI #3) | **accepted** §5.4 |
| G2 weight applied in the attribution prior (domain should) | **accepted** §5.4 |
| Adoption rule: both folds + paired bootstrap 90% LB > 0 (ext-PI #1, AI #N1) | **accepted** §8.2 |
| Greedy tiered tuning, arm handling first; cap = configurations evaluated (ext-PI #2, AI #N1) | **accepted** §5.5 |
| Minimum viable pipeline first; deferrable list (ext-PI #4) | **accepted** §12 |
| Seal TEST labels and raw DE outside the tree; commit LOCK; selected.yaml before final (AI #R1, #R2) | **accepted** §6 |
| G2 within-class and own-gene-stratified CMH reporting (AI #N4) | **accepted** for reporting. Tuning G2 on class (i) only is **rejected**: a regulator that is not expressed cannot be knocked down, so the pooled gain from G2 is real DepMap-derived knowledge, not an artefact. It is kept visible |
| ACH IDs asserted against Model.csv (AI #N4) | **accepted** §5.3 |
| NC2 restricted to DepMap-non-expressed regulators, evaluated without G2 (domain should) | **accepted** §7 |
| Hinge at copy-neutral (domain should) | **accepted** §5.5 |
| TP53 from 24Q4, rule Damaging > 0 or Hotspot > 0 (data #1) | **accepted** §4 |
| TP53-WT-only variant (domain should) | **rejected**: variant cap; TP53 is already a covariate |
| Regulator homozygous-deletion term; better ploidy imputation (domain could) | **deferred** |
| c₀ definition (ext-PI #6) | **accepted** (applied in v3.1 after round 3 caught the omission) |
| TP53 LOCO for chr17 regulators (ext-PI #7) | **modified**: keep TP53 as a covariate, because the mutation is upstream of 17p LOH, so it is a confounder rather than a mediator; dropping it for chr17 regulators is a sensitivity analysis. TP53 is pinned to 24Q4 |
| Ablations relative to the selection (ext-PI #9) | **accepted** §8.4 |

### Round 3 (fresh external PI; reviews/round3/external_pi_fresh.md). Dispositions:
| proposal | decision |
|---|---|
| #1 sign-rule contradiction | **accepted** §1, §6 |
| #2 c₀ | **accepted** §6 |
| #3 final score S_AB | **accepted** §5.6 |
| #4 resampling conventions; fixed top-K threshold; plain arm jackknife | **accepted** §6 |
| #5 dropped round-2 rows; TP53 pinned to 24Q4 | **accepted** §2, §13 |
| #6 terminology; #9 arm CN definition; #10 MVP call-set rule; #11 A3 notation; #12 CN-only regulators | **accepted** |
| #7 adoption-rule details | **accepted** §8.2 |
| #8 ablations; #13 TP53 chr17 | **accepted** (see round-2 rows) |
| #14 seeds via SeedSequence; #15 milestone duties; #17 operational exclusions | **accepted** |
| #16 TRAIN–TEST CN correlation diagnostic | **accepted** as a reported diagnostic |

The design is frozen after round 3. Further changes need a decision note in RESEARCH_TREE.md.
