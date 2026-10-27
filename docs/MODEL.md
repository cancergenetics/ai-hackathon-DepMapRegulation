# Model engine: implementation notes (DESIGN v3.1 §3–§5, §7)

Code: `src/depmaptr/{io,covariates,assoc,pipeline,gates,attribution,calibrate,baselines,simulate}.py`;
CLI: `scripts/run_model.py`; variants: `results/variants.yaml`; tests: `tests/test_model_*.py`, `tests/test_sim.py`.
This file lists the implementation choices where DESIGN leaves room. Anything not listed follows DESIGN literally.

## 1. Running

```
.venv/bin/python scripts/run_model.py --variant M1 --out data/results/scores/M1
.venv/bin/python scripts/run_model.py --variant focal+g1_t30_w0 --out data/results/scores/<run_id>   # composition
.venv/bin/python scripts/run_model.py --variant NC1 --perm-index 3 --out data/results/scores/NC1_3
.venv/bin/python scripts/run_model.py --variant B0 --seed 1 --out data/results/scores/B0_s1
.venv/bin/python scripts/run_model.py --variant M1 --sim 1 --out /tmp/x                              # SIM data
```
`--variant a+b` starts from `defaults` in `results/variants.yaml` and applies the `params` of a, then b.
Outputs follow docs/INTERFACES.md (`scores.npz`, `meta.json`) plus `manifest.json` (DESIGN §8.7).
Extra npz keys: `alpha` (A3 attribution weights), `t_gain` (hinge gain-slope t) and `call` (bool, MVP call
set, §7b; model runs only). A = B entries of `t` and `score` are NaN in every model run, NC1 included.
Regulators and targets are ordered by (chromosome 1–22, X; arm p < q; start).

## 2. Data and transforms (io.py)

- Panel order = `panel_s2.tsv` row order; CN and expression rows are taken from the matched parquets
  by ModelID (asserted complete). S1 = `in_S1`; `S2_noHCT116` drops ACH-000971.
- CN genes used for arms, segmentation and PCs: all CN columns on 1–22/X, not PAR/Y/MT, with an arm,
  finite in every panel line (18,998 on S2). Regulators = `variable_cn` and the same scope (12,141);
  targets = `measured` and the same scope (18,276).
- x = log2(clip(relCN, 2^-2, 2^3)); floor variant 2^-5. The homozygous-deletion indicator always uses raw
  relCN < 0.25 (≥ 3 lines; 691 targets on S2).
- **RINT variant (interpretation):** rank-based inverse normal Φ⁻¹((rank − 0.5)/n), average ranks for
  ties, applied per column to gene-level CN (x_A, x_B, cis) and to expression. **Arm CN and burden keep
  the log2 clip transform**, because the burden thresholds (|arm| > 0.3) are in log2 units; a literal
  "median of RINT values" would make them meaningless. Flagged for the PI.
- Context expression (G2, B4e) is read from the full matched expression matrix (ACH-000971 / ACH-001085);
  `io.context_model_check` asserts the IDs against Model.csv CellLineName (alphanumeric comparison, since
  Model.csv spells HEK 293 as "HEK-293").

## 3. Covariates Z (covariates.py)

- **Lineage:** for Lymphoid, Myeloid, CNS/Brain, Esophagus/Stomach, Soft Tissue, Bone, each primary
  disease with ≥ 15 lines becomes `Lineage::Disease`; the rest of the lineage stays as the parent level
  **even if < 15** (Myeloid 13, Bone 12, CNS/Brain 10 on S2). Other lineages with < 15 lines and missing
  lineage → `other` (47 lines). **32 levels** on S2 (31 dummies + the reference, the largest level, Lung).
- **Ploidy, WGD:** missing → median within the final lineage level (global median if none observed) +
  missingness indicator. `wgd_missing` is identical to `ploidy_missing` and is dropped as a duplicate.
- **MSI:** MSIScore > 20 (asserted equal to the panel's `msi_high`), missing → 0 + `msi_missing`.
- **TP53:** altered flag, missing → 0 + `tp53_missing` (identical to `msi_missing` on S2 → dropped).
- **Sex:** `Female` and `Unknown` (Unknown or missing) indicators, reference Male.
- **Media:** base-medium levels with ≥ 15 lines one-hot (reference RPMI); other levels and missing are
  pooled into `media[other]`. Serum-free flag, missing → 0 + `serum_free_missing`.
- **Growth pattern:** reference Adherent; Suspension, Mixed, Spheroid and `Unknown` (8 lines) indicators.
- **RNA technical:** `stranded`. **Platform dropped** (PI decision D3: 93% missing).
- Constant and exactly duplicated columns are dropped and logged in meta.json. On S2, Z has 56 columns.

## 4. Arms, burden and the design D per group

- arm_k = median of x over the CN genes on arm k (arm from the UCSC band at the gene midpoint, as in
  `gene_annotation.tsv`); arms with ≥ 20 genes (41 on S2). Xp/Xq are sex-centred (median within
  Male / Female / Unknown-or-missing). Autosomal arms agree with `depmap_arm_cn_derived.parquet` (tested).
- **Burden** = [fraction of autosomal arms with |arm| > 0.3, mean |arm|] over the autosomal arms that are
  in the arm block of the design:
  - LOCO (M1): group = chromosome c; arms ∉ c (incl. Xp/Xq for autosomal c); burden excludes c.
  - focal: one group; all 41 arms; burden over all autosomal arms (A's arm is adjusted anyway).
  - leave-one-arm-out: group = A's arm; all arms except it (A's other arm kept); burden excludes it.
    A regulator on an arm with < 20 genes (not an arm covariate) gets all arms.
  - CN-PCs: group = chromosome c; burden LOCO(c); top 20 PCs of the centred gene-level CN (active
    transform, X genes sex-centred) of chromosomes ≠ c, from the n×n Gram matrix (scores = U·√λ).
- D = [1, Z, burden, arm block]. Pivoted QR on column-normalised D (rank tolerance 1e-9 relative);
  cond of the retained R must be < 1e6 (AssertionError otherwise). S2 LOCO: rank 98 = all columns, cond ≈ 37.

## 5. Engine (assoc.py)

- Generic pair-specific-regressor engine: for each (A, B) the Gram matrix of [V_A, W_B, u_A, y_B] (all
  residualised on D) is built from GEMMs and broadcast vectors and swept elementwise. This is
  algebraically identical to DESIGN §5.1's "orthonormalise W̃_B into Q_B" formula; it also covers the hinge
  (V_A = gain part) and the co-expression term (u = y_g, V = x_g). float64 throughout, regulators chunked
  (512 rows).
- df = n − rank(D) − (#nuisance columns with a non-degenerate pivot) − 1. A nuisance pivot is degenerate
  if its residual sum of squares is ≤ 1e-10 × its value before the sweep; an inactive deletion indicator is
  an all-zero column and therefore never counted. So k_B and df are per target (and per pair if x_B is
  collinear with x_A, which only happens for same-chromosome pairs).
- Guards: NaN where √S_uu or √S_yy (after D, V, W) < 1e-3 × the centred raw norm, or df < 3. NaN → minimum
  score in evaluation (harness). On S2 M1 the NaNs are the 11,713 A = B pairs plus 3 more (same-chromosome).
- t = r √(df / (1 − r²)). p-values are not stored; `calibrate.z_from_t` uses log-space tails
  (with Wallace's approximation where scipy's log-tail underflows) for λ_GC.
- Known-answer test: statsmodels OLS on 240 random real pairs per arm handling (LOCO, focal; half with an
  active deletion indicator; chromosomes 1, 8, 17, X) and on toy data incl. a rank-deficient D, k = 3, hinge
  and co-expression forms; max relative error < 1e-8 (observed ~1e-14 on toy data).

## 6. Gates, attribution, final score

- **G1:** t_cis from y_A ~ x_A | D_group(A) (no deletion indicator; DESIGN §5.3), linear x even under the
  hinge, unpermuted under NC1. CN-only regulators (426 on S2) have t_cis = 0. s_A = −1 iff t_cis ≤ −3
  (3 regulators on S2).
- **G2:** expressed iff context log2(TPM+1) ≥ threshold; CN-only regulators are not expressed.
- **A3:**
  - Segmentation (`greedy_segments`) runs on x residualised on the **LOCO** D_c, whatever the arm-handling
    variant, so segments are a property of the CN data, not of the variant. Every in-scope CN gene of the
    chromosome is segmented and scored (t_gB for all of them), then α is normalised over all of them;
    only regulator rows are written out.
  - z_gB = t_gB (NaN → 0), with the pair's df in λ.
  - z_co,gB = atanh(r)·√(df − 1), with r the engine's partial correlation of y_g with y_B given (x_g, W_B, D)
    and that model's df. **CN-only genes (no expression column) get z_co = 0 with df_co = df − 1** (the df the
    co-expression model would have), so they contribute no co-expression evidence and never produce NaN
    (PI decision D5; before D5 they had df_co = NaN, which made α NaN for every gene of their segment).
  - π_g ∝ exp(t_cis,g / s)·w_G2(g) in log space (w_G2 only if G2 is active). If all π are 0, π is uniform.
  - `run` asserts that the final score has no NaN where t is finite.
- **Target standardisation:** per target, z-score of t across regulators on other chromosomes (NaN
  excluded, ddof = 1). Same-chromosome entries are transformed with the same mean/SD but are masked by
  the harness anyway.
- **S = t′ · w_G1 · w_G2 · α, always** (PI decision D5, amending DESIGN §5.6): w_G2 is applied as a multiplier
  for every regulator, including singleton segments and segments where all π are 0, and under A3 it
  additionally enters π. `t` in the score file is always the raw reduced-form t (the loss-slope t for the
  hinge).
- **Hinge:** x⁻ = min(x, 0), x⁺ = max(x, 0) at copy-neutral; for X regulators x is first sex-centred
  (median within sex group), because male X is at log2 ≈ −1 (review m-1). The linear part is unchanged
  because sex is in Z.

## 7. NC1 and baselines

- **NC1:** one within-lineage-level permutation per pass (strata = the 32 Z lineage levels), applied to
  the residualised regulator matrix (and V) of every group; Ỹ and W̃_B fixed. Seed: spawn index k of
  `SeedSequence(20261027)`.
- **B0:** uniform floats from `default_rng(SeedSequence(seed))`, seed default 20261027; the seed used is
  always written to meta.json `seed`.
- **B1:** Pearson t (df = n − 2) of x_A vs y_B. s_cis from the naive cis t.
- **B2:** Pearson t of y_A vs y_B; CN-only regulators NaN; s_cis = +1 everywhere (an mRNA statistic:
  lowering A predicts −sign(r)). cis_t holds the M1 cis t for reference only.
- **B3:** M1 cis t × SD(x_A) (signed; ranked by |·|). **B4:** SD(y_B). **B4e:** y_B in ACH-000971.
- **B5:** rank(|B3|) × rank(B4e) (average ranks), stored as float32 (products < 2.3e8; top-end spacing
  ≫ float32 resolution, so few artificial ties).
- **B6:** t of the arm CN of A's arm (LOCO design, same W_B); regulators on arms without an arm
  covariate (< 20 genes) are NaN.
- For B0, B3–B6 s_cis comes from the M1 cis t (target-agnostic baselines carry no sign information).

## 7b. MVP call set (DESIGN §12)

For model runs (not NC1, not baselines) `run_model.py` also runs `--nc1-perms` (default 1) NC1 passes of
the **same variant** in-process (indices 0..k−1; gates, A3 and target standardisation are switched off in
those passes, since they do not change t). Over trans pairs only, t* is the smallest |t| threshold
(1e-3 grid) with (#NC1 ≥ t* / n_perm) / (#real ≥ t*) ≤ 0.05, and `call = trans & |t| ≥ t*`
(`calls.py`). meta.json `calls` records t*, n_calls, the FDR estimate at t*, n_perm and the NC1 parameters.
Because the NC1 passes are always computed for the variant itself, the "fall back to M1's NC1" case never
arises. The call set is on t, not on the gated score.

## 8. Simulation (simulate.py; tests/test_sim.py)

600 lines; 6 autosomes × 2 arms × 60 genes + X (2 × 25); segments of 3–8 genes; 6 lineages; sex.
Planted: 12 focal regulators (one negative-cis, one on X) × 25 targets each, acting through A's mRNA, with
weak-cis co-segment neighbours. Confounders:
- 1q/3p co-selection in both directions with equal magnitude (so |arm|-based burden cannot absorb it), with
  a 3p arm-level driver of 20 targets;
- lineage L1 with 2q gains and an expression programme on 40 genes;
- a male-up autosomal programme (raw X CN tracks sex);
- homozygous deletions in 4p (expression collapses) co-occurring with amplification of a 5q partner segment;
- 6% CN-only genes, including a neighbour in half of the planted segments.

Focal losses are −0.9, so the relCN < 0.25 indicator marks only true homozygous deletions.

**Each confounder test has an ablation that must fail,** checked on seeds 1–3 in the suite, using monkeypatches
in `tests/test_sim.py::_ablate`:

| test | passes with | fails without |
|---|---|---|
| 1q → driver targets: FP < 2%; B1 > 10% | LOCO arm block | arm block removed: > 10% |
| 2q → lineage programme: FP < 2% | lineage in Z | lineage removed: > 30% |
| X regulators → sex programme: FP < 2% | sex in Z | sex removed: > 30% |
| 5q partner → deleted targets: FP < 5% | deletion indicator in W_B | indicator removed: > 30% |
| knockdown sign of the negative-cis regulator | s_A rule | s = +1: < 5% concordant |

Also tested:
- planted AUROC > 0.98, with B0 at chance;
- the planted X regulator's targets are recovered;
- A3 ranks the planted gene first in ≈ 90% of segments (chance ≈ 18%);
- A3 with w_co = 0.5 plus G2 adds no NaN, and gives 0 for every non-expressed regulator;
- the hinge on X matches statsmodels with sex-centred hinge regressors;
- λ_GC is 0.9–1.12 under the null and 0.85–1.12 for NC1;
- focal removes the arm-level LOCO signal.

Not simulated or tested: burden, TP53 and MSI confounding.

## 9. Real-data runs

All runs made before PI decision D5 were deleted. Only M1, NC1 and B0–B6 were regenerated with the fixed
code, on 2026-10-27. Tuning runs are driven by the PI through the harness.

Data: S2, n = 1,436; R = 12,141 regulators (426 CN-only); T = 18,276 targets; 221.9 M pairs, of which 209.9 M
are trans. Everything runs in a single process on Accelerate BLAS. In the table, "footprint" is the
`/usr/bin/time -l` peak memory footprint. macOS `ru_maxrss` (stored in meta.json) overstates it: 13.3 GB
for M1, 11.3 GB for NC1 and 5.4–6.6 GB for the baselines.

| run | compute s | end-to-end s | footprint GB | λ_GC (2 M entries) | NaN scores |
|---|---:|---:|---:|---:|---:|
| M1 (+ 1 in-process NC1 pass for calls) | 51.6 | 82 | 6.4 | 2.10 | 11,716 |
| NC1 (index 0) | 27.5 | 57 | 4.9 | 0.96 | 11,713 (A = B) |
| B0 | 2.7 | 32 | 3.0 | – | 0 |
| B1 | 4.1 | 34 | 3.1 | 4.75 | 0 |
| B2 | 4.8 | 34 | 3.0 | 16.7 | 7,822,128 (whole rows: 426 CN-only + 2 regulators with constant expression in S2) |
| B3 / B4 / B4e / B5 | 2.4–2.9 | 31–32 | 3.0–3.9 | – | 0 |
| B6 | 12.2 | 42 | 3.2 | 2.50 | 18,276 (1 regulator on an arm without an arm covariate) |

- M1 NaN = the 11,713 A = B pairs plus 3 guarded same-chromosome pairs.
- The B1 and B2 A = B entries are finite or infinite; they are same-chromosome entries, which the harness masks.
- M1 call set (§7b; 1 NC1 pass): t* = 3.016, with 10,673,076 trans pairs called (5.1% of trans pairs) at an
  estimated FDR of 0.0499.

**Calibration.** NC1 is slightly below 1, as residual permutation predicts: a permuted residual vector has
E[r²] ≈ 1/(n − 1), while t uses df = n − rank(D) − k_B − 1 ≈ 1,335, so E[t²] ≈ 1,335/1,435 ≈ 0.93.

M1's λ ≈ 2 is not an engine artefact:
- NC1 is calibrated, and λ ≈ 1 on the SIM null;
- the independent review (reviews/code/model_review.md) attributes most of the inflation to A's own-arm
  dosage, the DESIGN §1 LOCO caveat. In its check λ was 1.23 for focal and 2.49 for B6, and M1 rows correlated
  r ≈ 0.85 with B6 rows.

The Efron empirical null (DESIGN §5.2) is not implemented; it is deferrable per §12.
