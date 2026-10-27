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
Extra npz keys: `alpha` (A3 attribution weights) and `t_gain` (hinge gain-slope t).
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
  lineage → `other` (47 lines). 31 levels on S2, reference = largest (Lung).
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
    and that model's df. Genes without expression get z_co = 0.
  - π_g ∝ exp(t_cis,g / s)·w_G2(g) in log space (w_G2 only if G2 is active). If all π are 0, π is uniform.
- **Target standardisation:** per target, z-score of t across regulators on other chromosomes (NaN
  excluded, ddof = 1). Same-chromosome entries are transformed with the same mean/SD but are masked by
  the harness anyway.
- S = t′ · w_G1 · w_G2 · α (w_G2 omitted when A3 is active, where it enters through π). `t` in the
  score file is always the raw reduced-form t (the loss-slope t for the hinge).

## 7. NC1 and baselines

- **NC1:** one within-lineage-level permutation per pass (strata = the 31 Z lineage levels), applied to
  the residualised regulator matrix (and V) of every group; Ỹ and W̃_B fixed. Seed: spawn index k of
  `SeedSequence(20261027)`.
- **B0:** uniform floats from `default_rng(SeedSequence(seed))`, seed default 20261027.
- **B1:** Pearson t (df = n − 2) of x_A vs y_B. s_cis from the naive cis t.
- **B2:** Pearson t of y_A vs y_B; CN-only regulators NaN; s_cis = +1 everywhere (an mRNA statistic:
  lowering A predicts −sign(r)). cis_t holds the M1 cis t for reference only.
- **B3:** M1 cis t × SD(x_A) (signed; ranked by |·|). **B4:** SD(y_B). **B4e:** y_B in ACH-000971.
- **B5:** rank(|B3|) × rank(B4e) (average ranks), stored as float32 (products < 2.3e8; top-end spacing
  ≫ float32 resolution, so few artificial ties).
- **B6:** t of the arm CN of A's arm (LOCO design, same W_B); regulators on arms without an arm
  covariate (< 20 genes) are NaN.
- For B0, B3–B6 s_cis comes from the M1 cis t (target-agnostic baselines carry no sign information).

## 8. Simulation (simulate.py; tests/test_sim.py)

600 lines, 6 autosomes × 2 arms × 60 genes + X (2 × 25), segments of 3–8 genes, 6 lineages, sex.
Planted: 12 focal regulators (one negative-cis) × 25 targets each through A's mRNA; weak-cis co-segment
neighbours; 1q/3p arm co-selection with a 3p arm-level driver of 20 targets; lineage L1 with 2q gains and
an expression programme on 40 genes; homozygous deletions in 4p. Tested on 3 seeds: AUROC > 0.98;
1q → driver-target false positives under B1 (> 10% at |t| > 3) but not M1 (< 2%); lineage programme FPs
removed; sign rule correct incl. the negative-cis regulator; A3 ranks the planted gene first ~90% vs ~18%
chance; λ_GC within 0.9–1.12 under the null and 0.85–1.12 for NC1; focal removes arm-level LOCO signal.

## 9. Real-data run (S2, M1 + baselines, 2026-10-27)

n = 1,436; R = 12,141 regulators (426 CN-only); T = 18,276 targets; 221.9 M pairs (209.9 M trans).
Engine wall time: M1 25 s (one process, Accelerate BLAS; the whole process including loading, sha256 of
the inputs and writing the 1.8 GB npz takes 55 s). Baselines take 2–13 s of compute (~32–43 s end to end).
Peak RSS is ~5.5 GB for the baselines and ~11 GB for M1/NC1.

Calibration (2 M random entries; per-chromosome λ in meta.json):

| run | λ_GC | share with \|t\| > 5 |
|---|---:|---:|
| M1 | 2.10 | 0.27% |
| NC1 (index 0) | 0.96 (per-chromosome median 0.91) | 0.001% |
| B1 (naive) | 4.75 | 5.6% |
| B6 (arm only) | 2.50 | 0.57% |

The NC1 value slightly below 1 is expected for residual permutation: a permuted residual vector has
E[r²] ≈ 1/(n − 1), while t uses df = n − rank(D) − k_B − 1 ≈ 1,335, so E[t²] ≈ 1,335/1,435 ≈ 0.93.
M1's λ ≈ 2 is not an artefact of the engine: NC1 is calibrated, and on SIM the null gives λ ≈ 1. It reflects
genuine trans association structure that the covariates do not absorb, i.e. real effects plus
residual confounding. The empirical-null step of DESIGN §5.2 (Efron z′) is meant to handle this, and it
is not implemented here (it is deferrable per §12).

### Variant runtimes on S2 (scored only; nothing was evaluated against the benchmark)

| variant | engine s | end-to-end s | ru_maxrss GB | median per-chrom λ_GC |
|---|---:|---:|---:|---:|
| M1 (LOCO) | 25 | 55 | 11.3 | 1.96 |
| focal | 19 | 48 | 10.8 | 1.20 |
| loao | 34 | 64 | 10.6 | 1.90 |
| cnpcs20 | 29 | 58 | 11.6 | 2.00 |
| hinge | 61 | 92 | 11.8 | 1.45 (loss slope) |
| rint | 28 | 57 | 11.7 | 2.22 |
| floor5 | 25 | 55 | 11.4 | 1.94 |
| tstd | 28 | 57 | 11.7 | — |
| g1_t30_w0+g2_thr1_w0 | 26 | 55 | 11.4 | — (7,457 pass G1; 7,080 expressed in HCT116) |
| a3_r095_s2_w05+g2_thr1_w0 | 85 | 114 | 13.4 | — (5,126 segments at ρ = 0.95) |
| S1 (n = 999) | 20 | 50 | 10.5 | 1.63 |
| NC1 | 26 | 56 | 11.4 | 0.91 |

Every design was full rank with cond(D) ≤ 37. macOS `ru_maxrss` overstates the working set (the
`/usr/bin/time -l` peak footprint for M1 was 4.9 GB). In S1, 2 targets are constant among the 999 lines
and are NaN (guard). The focal variant's much lower λ is consistent with LOCO t being dominated by
arm-level dosage shared along A's arm (DESIGN §1 caveat).
