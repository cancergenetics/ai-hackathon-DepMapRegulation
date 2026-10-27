# DESIGN — DepMap CN → trans-expression effects (DepMapTR)

Version: v1 (2026-10-27, PI draft). Inputs: PROMPT.md (fixed goal), LITERATURE.md, DATA.md, docs/DATA_PROFILE.md.
Revision history and reviewer responses are in §12.

## 1. Goal and estimand

For a regulator gene A and a target gene B on different chromosomes, estimate the effect of A's copy
number on B's mRNA (β_AB). We want the total effect, both direct and indirect, of changing A's dosage in
a cancer cell, holding fixed the CN of B, B's arm, the rest of the genome's arm-level CN and
cell-intrinsic context. The tool outputs, for every trans pair (A, B), a signed effect, a test
statistic/p-value, and a ranking score. A pair is "called" at BH FDR ≤ 0.05.

A causal reading is only justified when A's CN changes A's own mRNA (the "instrument"). The CN signal is
shared by every gene in A's copy-number segment, so separating A from its co-segmented neighbours is a
central sub-problem (§5.3).

## 2. Data (details in DATA.md / docs/DATA_PROFILE.md)

- DepMap Public 26Q3 (inferred release; checksums in DATA.md): `OmicsCNGeneWGS` (linear relative CN,
  ploidy-normalised) and `OmicsExpressionTPMLogp1HumanProteinCodingGenes` (log2(TPM+1)). Use the
  `IsDefaultEntryForModel` profile per model. 1,724 models are matched (1,509 cell lines, 215 organoids).
- Covariates: Model.csv (lineage, sex, GrowthPattern), ModelCondition (media), DepMap 24Q4 PureCN
  (ploidy, WGD, CIN, LoH fraction, aneuploidy; 1,037 matched lines), MSIScore (1,311), and arm-level CN
  derived from gene CN (all models).
- Benchmark: PIE X-Atlas/Orion `de/*.parquet` (HF commit abfa3da), HCT116 and HEK293T. These are
  pre-computed perturbation–mRNA associations only (Mann–Whitney p, BH within target, LFC).

## 3. Sample set and preprocessing

- **Primary sample set S1:** cell lines (not organoids) with complete PureCN covariates (n ≈ 1,037;
  this matches the brief's "~1,000"). Variant S2 adds all 1,509 cell lines, with missingness indicators
  and imputed ploidy/WGD/CIN. The S1 vs S2 choice is a tunable compared on training regulators only.
  Organoids are excluded throughout: different culture system, no covariates.
- **CN:** x = log2(max(relCN, 0.05)), winsorised at the 0.5/99.5 percentiles per gene. Variant: rank-based
  inverse-normal transform.
- **Expression:** log2(TPM+1) as provided. A DepMap target is "measured" if present and SD > 0 in S1.
- **Regulators analysed by the tool:** all autosomal genes with variable CN, defined a priori as ≥ 20
  lines in the sample set with |log2 relCN| ≥ 1. X/Y regulators are excluded, because their CN is
  dominated by sex and X-inactivation. X/Y targets are kept, with sex as a covariate.

## 4. Covariates (justification)

Shared covariates Z enter every model:
| covariate | why |
|---|---|
| lineage (OncotreeLineage, levels < 10 merged to "other") | strongest expression confounder; lineage-specific arm events (e.g. CRC 18q−/20q+) create cross-chromosome CN–expression correlations |
| ploidy, WGD | ploidy-normalised CN hides absolute dosage; WGD alters global transcription and cell size |
| CIN / aneuploidy score / LoH fraction | the generic aneuploidy-stress axis (proliferation, ISR, p53) would otherwise make many CN-variable regulators "regulate" the same stress genes |
| MSI score (MSI-high flag) | MSI lines are near-diploid, with a distinct expression programme; HCT116 itself is MSI-H |
| sex | X/Y target dosage |
| media (base medium + serum-free flag, levels ≥ 10 lines) and growth pattern (adherent/suspension/mixed) | culture-driven expression programmes |
| purity | not available; cell lines ≈ 1 (assumption A5 in DATA_PROFILE) |

Pair-specific covariates: **CN_B** (the target's cis CN; required by the brief) and **arm-level CN of all
arms except those on A's chromosome** ("leave-one-chromosome-out", LOCO). LOCO arms absorb correlated
arm events on other chromosomes, which the brief warns about. A's own arm is excluded, because it is
collinear with CN_A and would remove the signal we want. B's own arm is included in the LOCO set,
covering B's regional dosage.

Deliberately **not** used: expression PCs/PEER factors. They absorb genuine CN-driven trans hubs
(Fehrmann 2015; Kang 2008; Dahl 2019 in LITERATURE.md). A CN-derived-PC variant is allowed. A
proliferation score derived from expression is reported only as a sensitivity analysis.

## 5. Statistical model

### 5.1 Core model (M1: covariate-adjusted reduced form)
For each pair:  y_B = β_AB·x_A + γ·x_B + Σ_k δ_k·arm_k(k ∉ chr(A)) + Z·η + ε.
β_AB is estimated by Frisch–Waugh–Lovell. Within each regulator chromosome c (22 groups):
1. Residualise Y (all targets) and X (CN of all genes) on [Z, arms∉c] (one QR per chromosome).
2. Compute the partial correlation of x_A and y_B given x_B in closed form, from three residual
   correlation matrices: r(x_A, y_B), r(x_A, x_B) and r(x_B, y_B). The last is a per-target vector.
3. Convert to t with df = n − p − 2. p-values come from t, then BH over all trans pairs.

This conditions exactly on Z ∪ arms∉c ∪ {x_B}. Cost is O(n·G_A·G_B) per chromosome via BLAS. That is
about 10^8–10^9 pair tests in minutes, chunked to stay inside 32 GB.

### 5.2 Instrument gating and IV scaling (M2)
- The cis effect of A is the partial regression of y_A on x_A given Z and arms∉chr(A). Regulators with
  no positive cis effect (cis t < τ_cis) cannot transmit CN to B through A's own product. They are
  down-weighted or excluded, with τ_cis tuned.
- The IV (Wald / 2SLS) estimate is β_IV = β_AB / β_cis,A: the effect per unit of A mRNA. It is reported
  for interpretation. Ranking keeps the reduced-form t, because the t-statistics are identical.

### 5.3 Co-segment attribution (M3)
Every gene in A's segment shares the CN signal, so the reduced form cannot tell A from its neighbours.
Two DepMap-internal ways to break the tie, both tuned on training regulators:
- (a) **Mediation evidence:** the partial correlation of y_B with A's non-CN expression variation
  (y_A residualised on x_A, Z, arms), given x_B. If A causes B, A's expression should predict B beyond
  shared CN. The combined score is the reduced-form evidence × mediation evidence (e.g. sum of signed
  z's when the signs agree).
- (b) **Segment-level discovery plus within-segment attribution:** collapse highly correlated adjacent
  genes (|r(CN)| > 0.95) into segments, test segment CN vs B, and attribute to genes by cis strength
  and mediation evidence (a SuSiE-like posterior over genes in the segment; simpler softmax variant
  first).

### 5.4 Robustness variants (tuned / compared on training regulators)
- Robust: rank/inverse-normal transformed Y and X (Spearman-type partial correlations).
- GLS / LMM: whiten by a CN-similarity kinship matrix (one eigendecomposition) for residual population
  structure.
- Lineage-stratified meta-analysis (within-lineage estimates combined by inverse variance) vs pooled
  with lineage fixed effects.
- Context weighting toward HCT116-like lines (colorectal / MSI-H). A secondary variant only, since it
  narrows generality.
- Empirical null: calibration checked by within-lineage permutation of CN profiles (whole CN vector
  per line, preserving CN structure).

## 6. Benchmark definitions (fixed before any model is run)

- **Cell lines:** HCT116 is primary (it is a DepMap model and a colorectal MSI-H carcinoma; it has 2×
  the targets and a plausible null level among non-expressed targets). HEK293T is secondary: not a
  DepMap model, and 83% of its trans hits come from non-expressed targets (DATA_PROFILE §11).
- **Positive (direction-agnostic):** trans pair with PIE `p_adj < 0.05` (BH within target, as released).
  Sensitivity analysis: `p_adj < 0.05 & |LFC| ≥ 0.25` (LFC winsorised to ±5).
- **Sign:** sign(log2FC) of the benchmark; sign(β_AB) of the prediction. **Concordant** iff
  sign(β_AB) = −sign(LFC_B|kd A): CRISPRi lowers A, so a positive CN→expression effect implies B goes
  down on knockdown. This is unit-tested on toy data.
- **Universe (primary, literal to the brief):** trans pairs (chrom(A) ≠ chrom(B), both placed) where:
  - A is a target in that cell line's DE table, maps 1:1 to DepMap, is autosomal, and has variable CN
    (as in §3, computed in S1);
  - B is consistently tested in that line (control mean ≥ 0.05, i.e. tested for every target) and
    measured in DepMap (§3);
  - A ≠ B.
- **Universe (secondary):** the primary universe restricted to regulators whose own gene is expressed
  and knocked down (own feature tested, p_adj < 0.05, LFC < 0). Non-expressed targets form an
  empirical-null regulator set: a good method should show OR ≈ 1 on them.
- **Held-out split (frozen before any tuning; seed 20261027):** split genomic regions, not genes, so CN
  neighbours cannot leak between splits. All autosomal genes are binned into 5 Mb bins; one third of
  bins are assigned to TEST within each chromosome (stratified by chromosome). A regulator's split is
  its bin's split. The split file is committed (`results/split/regulator_split.tsv`) before any model
  is evaluated. TEST labels are not read by any tuning code (enforced in code: the evaluation function
  refuses TEST unless `--final`).
- **Metrics:**
  - Primary: Fisher exact OR and one-sided p for the top K = 0.1% of universe pairs ranked by the
    method's score, HCT116 TEST universe. Ranking-based, so methods are compared at equal call
    volume, independent of calibration.
  - Secondary: OR at K ∈ {0.01%, 1%} and at the method's BH-5% call set; sign-concordant OR (positives
    redefined as benchmark-positive and concordant); concordance rate among overlapping pairs, with a
    binomial test vs 0.5; regulator-stratified Cochran–Mantel–Haenszel OR (guards against
    hub-dominance); per-regulator AUROC averaged over regulators with ≥ 10 positives; bootstrap 95% CIs
    resampling 5 Mb bins.
  - Runtime (wall clock, CPU threads) is reported for every method.

## 7. Baselines and controls

| id | description | purpose |
|---|---|---|
| B0 | random scores (seeds list) | null: OR ≈ 1 |
| B1 | naive Pearson r(x_A, y_B), no covariates | quantifies covariate value |
| B2 | co-expression r(y_A, y_B), no CN | is CN adding anything over co-expression? |
| B3 | regulator-hub-only score: rank pairs by a DepMap-only regulator property (cis strength × CN variance), target-agnostic | degree-bias baseline |
| B4 | target-hub-only: target expression variance across DepMap | degree-bias baseline |
| NC1 | M1 with CN profiles permuted across lines within lineage | negative control: OR ≈ 1 |
| NC2 | the secondary-universe non-expressed regulators | biology null |
| PC1 | cis pairs (A = B, own CN→own mRNA) | positive control for the CN pipeline (not part of metrics) |
| PC2 | HEK293T hits predicting HCT116 hits for shared regulators | empirical ceiling / reference OR |
| SIM | toy simulation with planted trans effects, arm confounders, lineage confounders and co-segment neighbours | verifies recovery, calibration, sign convention and that confounders are removed |

## 8. Rules guarding against overfitting and AI-agent failure modes

1. Benchmark labels are used only to choose among a small, pre-listed set of hyperparameters and
   variants (§5) on TRAIN regulators. They are never used as features (e.g. no learned target-hub
   prior from Perturb-seq hits). DepMap CRISPR essentiality, LINCS, other Perturb-seq and
   perturbation-trained embeddings are excluded from the primary method.
2. Every tuning run is logged to `results/tuning_log.tsv` (variant, params, TRAIN metrics, runtime,
   git hash). The number of variants tried is reported in REPORT.md.
3. TEST is evaluated once per final method (plus pre-declared baselines). Any re-run must be logged
   and justified.
4. Metrics, labels, universe and positives are frozen in `src/depmaptr/benchmark.py` with tests.
   Changes after freezing need a decision note and must be reported.
5. Each result relied on is reviewed by an agent that did not produce it.
6. Seeds list: [20261027, 1, 2, 3, 4] (split uses 20261027 only).

## 9. Software layout

`src/depmaptr/`: `io.py` (load matched matrices + covariates), `covariates.py`, `assoc.py`
(chunked FWL partial correlations), `attribution.py`, `benchmark.py` (universe, labels, metrics,
split guard), `simulate.py`. `scripts/`: pipeline entry points with CLI args. `tests/`: unit tests
plus simulation tests. `results/`: committed small outputs (metrics tables, logs, plots). The big
pair-score matrices go in `data/results/` (ignored).

## 10. Compute plan

1,037 × ~18.5k matrices. Scoring all ~9k variable-CN regulators × ~15k targets (1.4×10^8 pairs) is a
few BLAS GEMMs per chromosome: minutes on 10 cores. Store top pairs per regulator (sparse) plus full
float16 scores for universe regulators only.

## 11. Milestones

M3 harness + split + toy tests → M4 baselines + M1 → M5 variants (M2, M3, robustness) on TRAIN → M6
final TEST evaluation → M7 REPORT.md.

## 12. Review log
(Filled in after each review round.)
