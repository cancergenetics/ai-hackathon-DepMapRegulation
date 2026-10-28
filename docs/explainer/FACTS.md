# DepMapTR — verified fact sheet for explainer writers

Every number below was re-read from the file named in the "source" column on 2026-10-28 (or recomputed, where
marked "recomputed"). "Split" says whether the number comes from the TRAIN regulators (used for tuning), the
held-out TEST regulators (one locked run, `results/final/`, commits 603a7f5 → e9085c7), or is label-free
(no benchmark labels involved). Rounding follows REPORT.md unless stated.

Abbreviations: OR = odds ratio of Perturb-seq hits among the top 0.1 % of ranked pairs vs the rest;
perm p = regulator-label permutation p (1,000 permutations on TEST); jk CI = delete-one-arm jackknife 95 % CI;
A = regulator (the gene whose copy number varies / is knocked down); B = target (the gene whose mRNA responds).

---

## 1. Data sizes (label-free)

| fact | value | source |
|---|---|---|
| DepMap release | Public 26Q3, **inferred** (portal md5s not verifiable) | DATA.md; docs/DATA_PROFILE.md §1 |
| matched CN+RNA models | 1,724 (1,509 cell lines, 215 organoids) | docs/DATA_PROFILE.md §2 |
| primary sample set S2 (cell lines; no engineered models; one per patient) | **1,436** | `data/results/scores/g1_t100_w0__a3_r095_s2_w0/meta.json` `n_samples`; `data/processed/panel_s2.tsv` |
| S1 (PureCN-complete subset of S2; sensitivity only) | 999 (drops HCT116 itself) | docs/DATA_PROFILE.md §12; runtime_table `S1` n_samples 999 |
| PureCN ploidy/WGD missing in S2 (imputed + indicator) | 437 / 1,436 (30 %) | docs/DATA_PROFILE.md §12 |
| in-scope CN genes (autosomes + X, finite in all 1,436 lines; used for arms and segments) | 18,998 | docs/MODEL.md §2 |
| regulators scored (frozen `variable_cn`: ≥ 20 S2 lines with \|log2 relCN\| ≥ 1, X sex-centred) | **12,141** (11,425 autosomal + 716 X) | meta.json `n_regulators`; DATA_PROFILE §12 |
| … of which CN-only (no expression column) | 426 | meta.json `n_regulators_without_expression` |
| targets scored (frozen `measured`) | **18,276** | meta.json `n_targets` |
| all pairs / trans pairs (different chromosomes) | 221,888,916 / **209,943,680** | meta.json `n_pairs`, `n_trans_pairs` |
| shared covariate matrix Z | 56 columns (32 lineage levels → 31 dummies; ploidy, ploidy_missing, wgd, msi_high, msi_missing, tp53_altered, sex ×2, serum_free ×2, stranded, 10 media dummies, 4 growth-pattern dummies) | meta.json `diagnostics.Z_columns` |
| arm-level CN covariates | 41 arms (≥ 20 genes each); 39 used per chromosome group for chr1 (arms of A's chromosome excluded) | meta.json `diagnostics.n_arms`; recomputed (worked_example.json `design_columns` 98) |
| design D_c per regulator chromosome | 98 columns, full rank, cond ≈ 37 (chr1) | meta.json `diagnostics.groups.1` |
| targets with an active homozygous-deletion indicator (relCN < 0.25 in ≥ 3 lines) | 691 | meta.json `diagnostics.n_targets_with_homdel_indicator` |
| genome-wide inflation λ_GC of M1 t (NC1 null ≈ 0.96) | 2.10 | docs/MODEL.md §9 |
| A3 co-segments (ρ = 0.95), all chromosomes | 5,126 (chr1: 483) | meta.json `diagnostics.n_segments_by_chrom` (summed) |
| regulators passing G1 (t_cis ≥ 10) | **5,403 / 12,141** | meta.json `diagnostics.n_g1_pass`; recomputed from `scores.npz` `cis_t` |
| regulators with s_A = −1 (t_cis ≤ −3) | 3 | recomputed from `scores.npz` `s_cis` |
| regulators with t_cis = 0 | 428 (426 CN-only + 2 with constant expression) | recomputed from `scores.npz` `cis_t` |
| trans score entries exactly 0 in the selected method (G1-failed rows) | 116,452,204 of 209,943,680 (55 %) | recomputed, `worked_example.json` step7 |
| MVP call set (\|t\| ≥ t* = 3.016, NC1-FDR 5 %, 1 NC1 pass) | 10,673,076 trans pairs (5.1 %), FDR 0.0499 | meta.json `calls` |

## 2. Benchmark (PIE X-Atlas/Orion CRISPRi Perturb-seq, pre-computed DE tables, HF commit abfa3da)

Terminology trap: in PIE, `target` = the **perturbed** gene (our regulator A) and `feature` = the measured gene
(our target B).

| fact | HCT116 (primary) | HEK293T (secondary) | split | source |
|---|---:|---:|---|---|
| perturbations in DE table | 8,987 | 4,297 | — | `data/processed/bench/manifest.json` `regulator_counts.pie_targets` |
| universe regulators (mapped, variable CN, autosome/X) | **5,623** | 2,688 | both | manifest `universe_regulators` |
| universe targets (control mean ≥ 0.05, DepMap measured) | **10,209** | 11,000 | both | manifest `universe_targets` |
| universe trans pairs | 54,357,451 | 27,984,730 | both | manifest `n_pairs` |
| primary positives (p_adj < 0.05, BH within perturbation) | 115,601 | 52,510 | both | docs/BENCHMARK.md §1 |
| perturbations whose own gene is expressed/tested | 3,838 (43 %) | 436 (10 %) | — | DATA_PROFILE §10 |
| share of all trans hits from perturbations of **non-expressed** genes | 12.7 % | **83.0 %** | — | DATA_PROFILE §10–§11 |
| same, recomputed on TRAIN (own-class iii regulators' share of TRAIN positives) | 10.0 % | 87.8 % | TRAIN | recomputed from `positives_*_TRAIN.parquet` + `universe_regulators_*.tsv` `own_class` |
| hit concentration: share of trans hits from the top 1 % / 10 % of perturbations | 55 % / 92 % | 48 % / 94 % | — | DATA_PROFILE §10 |
| same on HCT116 TRAIN regulators | 56 % / 96 %; 52.5 % of TRAIN regulators have 0 hits; max 2,285 (POLR1G) | — | TRAIN | recomputed from `positives_hct116_TRAIN.parquet` |
| on-target LFC median (TRAIN) | −27.6 | −26.7 | TRAIN | manifest `median_on_target_lfc_train` (LFC magnitude is a pseudocount artefact; only sign is used) |
| HEK293T is not a DepMap line (HEK293 ACH-001085 is used as its context) | | ✓ | — | DATA_PROFILE §2, §11 |

## 3. Held-out split (frozen before tuning; seed 20261027)

| fact | value | source |
|---|---|---|
| unit | whole chromosome arm; 42 arm units, sorted by HCT116 regulator count, one random TEST arm per consecutive triplet | docs/BENCHMARK.md §2; `results/split/arm_units.tsv` |
| TEST arms (14) | 17q, 7q, 14q, 6q, 3p, Xp, 7p, 22q, 19p, 4q, 5p, 12q, 18p, 21q | BENCHMARK.md §2 |
| TRAIN arms | 28 (two folds of 14 for the adoption rule) | `results/split/train_folds.tsv` |
| HCT116 TRAIN: regulators / pairs / K (0.1 %) / positives | 3,806 / 36,664,698 / 36,665 / 81,258 (own class i/ii/iii 1,388/171/2,247) | BENCHMARK.md §2; `results/eval/g1_t100_w0__a3_r095_s2_w0/hct116_TRAIN.json` |
| HCT116 TEST: regulators / pairs / K / positives | 1,817 (32 %) / 17,692,753 / 17,693 / 34,343 | `results/final/g1_t100_w0__a3_r095_s2_w0/hct116_TEST.json` `metrics.n_*` |
| HEK293T TRAIN / TEST regulators | 1,890 / 798 | BENCHMARK.md §2 |
| sealing | TEST labels and raw DE tables moved to `~/DepMapTR_sealed`, chmod a-rw, sha256 in `results/FROZEN.sha256`; one run, LOCK committed before any TEST file was read | BENCHMARK.md §5; `results/final/LOCK` (started 23:08:15, finished 23:12:46, status COMPLETED) |
| split leakage | 10 chromosomes have one arm in each split; median max \|r\| of a TEST regulator's CN with any TRAIN regulator 0.19; Xp–Xq arm r = 0.55 | BENCHMARK.md §9 |

## 4. The final (selected) method — checked against the code

Run id `g1_t100_w0__a3_r095_s2_w0` (`results/final/selected.yaml` `selected`), params in its `meta.json`.

| component | exact definition in code | file:function |
|---|---|---|
| CN transform | x = log2(clip(relCN, 2⁻², 2³)) (relCN = DepMap linear relative CN, ploidy-normalised) | `io.py:cn_transform` |
| shared covariates Z | lineage (32 levels; heterogeneous lineages split by disease ≥ 15 lines), ploidy & WGD (lineage-median imputed + 1 missing indicator), MSI-high (MSIScore > 20) + missing, TP53 altered, sex (Female, Unknown), base medium (levels ≥ 15, ref RPMI), serum-free + missing, growth pattern, RNA `stranded`; constant/duplicate columns dropped (wgd_missing, tp53_missing) → 56 columns. Purity not available. | `covariates.py:build_Z`, `lineage_levels` |
| arm CN | median of x over the arm's genes (arms ≥ 20 genes, 41 arms); Xp/Xq sex-centred | `covariates.py:arm_cn` |
| LOCO design D_c for regulator chromosome c | [1, Z, burden(c), arm CN of every arm **not** on c]; burden = [fraction of autosomal arms not on c with \|arm\| > 0.3, mean \|arm\|] | `covariates.py:DesignBuilder.design`, `burden` |
| M1 reduced-form association | OLS t of x_A in y_B ~ D_c + x_B + 1(relCN_B < 0.25, only if ≥ 3 lines) + x_A, computed by Frisch–Waugh–Lovell (pivoted-QR residualisation on D_c, then a batched sweep over the cis basis); df = n − rank(D_c) − #active cis columns − 1 (= 1,336 typically) | `pipeline.py:group_t` → `assoc.py:chunked_pair_partial`/`pair_partial`, `Residualiser` |
| G1 cis-instrument gate | t_cis = partial t of y_A on x_A given D_c (no deletion indicator); **pass iff t_cis > 0 and t_cis² ≥ 100** (i.e. t_cis ≥ 10, first-stage F ≥ 100); w_G1 = 1 if pass else **0**; CN-only regulators t_cis = 0 → fail | `pipeline.py:run` (diag_partial), `gates.py:g1_pass`, `g1_weight` |
| knockdown sign rule | s_A = −1 iff t_cis ≤ −3, else +1 | `gates.py:s_cis` |
| A3 segmentation | per chromosome, CN of **all** 18,998 in-scope CN genes (not only regulators) residualised on the LOCO D_c; greedy left-to-right: a gene joins the current segment iff corr(gene, running segment mean) ≥ ρ = 0.95 | `attribution.py:greedy_segments` (called in `pipeline.py:run`) |
| A3 attribution (SuSiE L = 1) | for each target B, over the m genes g of a segment: α_{g\|B} ∝ π_g · ABF(t_gB); log ABF(z) = ½ log(1 − λ) + ½ λ z², λ = σ₀²/(σ₀² + 1/df), σ₀ = 0.05 fixed; π_g ∝ exp(t_cis,g / s), s = 2; co-expression weight w_co = 0 (off); singletons α = 1; computed in log space | `attribution.py:segment_alpha`, `log_abf` |
| final score | **S_AB = t_AB · w_G1(A) · α_{A\|B}** (no G2, no target standardisation); A = B and same-chromosome entries are not ranked | `pipeline.py:compose_score` |
| ranking | \|S\| over trans pairs; NaN → minimum; ties at the K boundary broken by a seeded hash | `benchmark.py:Evaluation.topk` |
| predicted knockdown direction of B | **−sign(S_AB) · s_A** (0 if S = 0 or NaN; 0 never counts as concordant) | `benchmark.py:Evaluation.pred_kd_sign`; `gates.py:predicted_kd_sign` |
| call set (separate from ranking) | trans pairs with \|t\| ≥ t*, t* the smallest threshold with NC1-permutation FDR ≤ 5 % (1 in-process NC1 pass of the same variant, lineage-stratified residual permutation) | `scripts/run_model.py`, `calls.py:t_star`, `call_matrix`, `pipeline.py:nc1_permutation` |
| context layer (reported, **not** selected) | S = sign(t_AB) · ρ_P · ρ_T(B)²; ρ_P = average-rank percentile of \|S_selected\| among all finite trans entries of the 12,141 × 18,276 matrix (zeros included; 0 where S = 0); ρ_T = percentile of B's DepMap HCT116 (ACH-000971) log2(TPM+1) among the 18,276 targets | `pipeline.py:context_fusion`, `trans_percentile`, `target_context_percentile` |
| engine correctness | matches statsmodels OLS to < 1e-8 on 480 real pairs (tests); worked example: \|Δt\| < 1e-14 for 3 refitted pairs; recomputed S row equals committed `scores.npz` to 5e-7 | docs/MODEL.md §5; `worked_example.json` |

## 5. Main TEST results (HCT116 unless stated; `results/final/summary.md` / `hct116_TEST.json`)

### 5.1 Pooled OR at top 0.1 % (K = 17,693 of 17,692,753 pairs)

| method | OR | jk 95 % CI | perm p | degree-perm p | OR@0.01 % / 1 % | score time s |
|---|---:|---|---:|---:|---|---:|
| **selected CN-only (M1 + G1 + A3)** | **3.85** (131 hits in top-K) | 1.33–11.2 | **0.068** | 0.061 | 2.92 / 4.77 | 67.2 |
| ablation: minus G1 (M1 + A3) | 1.46 | 0.91–2.33 | 0.229 | 0.206 | 0.87 / 2.02 | 66.1 |
| ablation: focal arm handling | 4.50 | 1.19–17.1 | 0.083 | 0.018 | 2.92 / 5.49 | 51.5 |
| context layer (selected + I1 γ = 2) | **18.4** | 10.9–30.9 | **0.015** | 0.021 | 16.5 / 13.4 | 82.7 |
| B0 random | 1.17 | 0.96–1.41 | 0.145 | 0.171 | | 2.7 |
| B1 naive Pearson x_A vs y_B | 1.14 | 0.31–4.12 | 0.723 | | | 4.1 |
| B2 co-expression y_A vs y_B | 6.03 | 2.71–13.4 | 0.095 | 0.332 | | 4.8 |
| B3 regulator hub (cis t × CN SD; target-agnostic) | 0.029 | (degenerate) | 0.666 | | 0.145 / 10.6 | 2.4 |
| B4 target expression SD | 5.73 | 4.23–7.76 | 0.528 | | | 2.7 |
| B4e target expression in DepMap HCT116 | 7.39 | 5.55–9.83 | 0.531 | | | 2.9 |
| B5 hub prior (rank B3 × rank B4e) | **22.8** | 3.27–159 | 0.043 | 0.035 | | 2.4 |
| B6 arm-only (A's arm CN) | 0.611 | 0.35–1.07 | 0.989 | | | 12.2 |
| NC1 permuted CN (negative control) | 1.28 | 0.48–3.40 | 0.163 | | | 27.5 |
| C2 context fusion on NC1 t | 10.6 | 8.57–13.1 | 0.148 | | | 42.8 |
| PC2 HEK293T hits → HCT116 hits (598 shared TEST regulators; empirical ceiling) | **181** | 75.2–435 | 0.001 | | | — |

Fisher p for the selected method is 2.0 × 10⁻³⁶ but assumes independent pairs (rejected by design).

### 5.2 Paired comparisons (Δlog-OR, paired jackknife 95 % CI; `summary.md` "paired")

| comparison | Δlog-OR | CI | |
|---|---:|---|---|
| selected vs minus-G1 (value of G1 on TEST) | **+0.97** | 0.30–1.64 | significant |
| selected vs B0 | +1.19 | 0.13–2.26 | significant |
| selected vs NC1 | +1.10 | −0.66–2.85 | **not** significant |
| selected vs B2 / B4 / B4e / B5 | −0.45 / −0.40 / −0.65 / −1.78 | all CIs include 0 | |
| selected vs focal | −0.16 | −0.69–0.37 | |
| context layer vs B5 | −0.22 | −2.29–1.86 | not different |
| context layer vs C2 | **+0.55** | 0.07–1.03 | significant |
| context layer vs selected | +1.56 | 0.71–2.42 | |
| leave-HCT116-out vs selected | −0.016 | −0.08–0.05 | (REPORT: "−0.02") |

### 5.3 Within-regulator information (CMH, AUROC, C1)

| metric | value | source field |
|---|---|---|
| CMH OR, strata regulator × target pan-DepMap expression decile | **0.94** (0.77–1.15), p 0.54 | `cmh_regulator_x_exprdecile` |
| CMH, strata regulator × HCT116 expression decile (D6 diagnostic) | 0.99 (0.81–1.21) | `cmh_regulator_x_context_exprdecile_D6` |
| CMH, strata own-class × expression decile | 1.80 (1.51–2.14) | `cmh_ownclass_x_exprdecile` |
| mean per-regulator AUROC (179 regulators with ≥ 10 hits) | **0.500** | `auroc` |
| C1 pair shuffle (within regulator, 50 HCT116-expression bins) | 3.76 vs 3.85; Δlog-OR +0.023 (−0.75 to 0.79) | summary.md "C1" (mode "score (approximation)", equivalent for a CN-only score) |
| C1 for context layer (exact mode) | 17.3 vs 18.4; Δlog-OR +0.06 (−0.06 to 0.19) | summary.md "C1" |
| share of top-K pairs from the 10 busiest regulators | 26 % | `top10_regulator_share` |
| hub-excluded (drop top 1 % regulators [19] and targets [103] by degree, re-rank) | 2.69, Fisher p 7 × 10⁻⁶ (context layer 14.8; B5 11.8) | `hub_excluded` |

### 5.4 Sign concordance (CN loss ↔ knockdown), selected method

| fact | value | source field |
|---|---|---|
| overlapping pairs (top-K ∩ hits) | 131 | `sign.n_overlap` |
| concordance | **0.718** | `sign.concordance` |
| null 1: pooled independence c₀ = qu + (1−q)(1−u) (q = 0.588 predicted down, u = 0.733 observed down) | 0.541 | `sign.c0` |
| null 2: within-regulator LFC-sign permutation, mean | **0.647** | `sign.sign_perm_null_mean` |
| permutation p (vs null 2) | **0.012** | `sign.sign_perm_p` |
| sign 2×2 OR / sign-concordant OR | 7.02 / 5.84 | `sign.sign_or`, `sign_concordant_or` |
| regulators with ≥ 5 overlaps / regulator-weighted concordance | 6 / 0.654 | `n_regulators_weighted`, `regulator_weighted_concordance` |
| minus-G1 ablation | 0.72 vs c₀ 0.49, n = 50, p 0.006 | summary.md "sign" |
| context layer | 0.549 vs 0.507, n = 601, p 0.22 | summary.md "sign" |
| NC1 | 0.364 vs 0.417, n = 44, p 0.52 | summary.md "sign" |
| robust in: noHCT116 (p 0.004), CN-outlier (0.007), hidden factors (0.012), TP53-chr17 (0.005), paralog exclusion (0.004), class i (0.008) | | summary.md "sign" |
| **not** robust in: S1 (0.586 vs 0.467, p 0.13), no-split-sibling (0.50 vs 0.48, n = 32, p 0.48) | | summary.md "sign" |

### 5.5 Sensitivity analyses and subsets of the selected method (HCT116 TEST OR@0.1 %)

| analysis | OR | jk CI | perm p |
|---|---:|---|---:|
| S1 sample set (n = 999) | 3.40 | 1.64–7.06 | 0.105 |
| leave-HCT116-out of DepMap | 3.79 | 1.38–10.4 | 0.069 |
| regulator-conditional hidden factors (k = 10) | 2.90 | 1.05–7.99 | 0.188 |
| CN-outlier lines excluded (15 lines) | 3.49 | 1.38–8.85 | 0.097 |
| TP53 covariate dropped for chr17 | 3.73 | 1.36–10.2 | 0.071 |
| paralog pairs excluded | 3.82 | 1.34–10.9 | 0.069 |
| positives p_adj < 0.01 / \|Ueffect\| ≥ 0.2 | 4.05 / 2.65 | | |
| own-class (i) regulators (own gene knocked down; 723) | 1.61 | 0.78–3.34 | 0.111 |
| own-class (iii) regulators (own gene not expressed; 1,030) | 0.23 (1 hit) | | 0.516 |
| NC2 (class iii and not expressed in DepMap HCT116): selected / ungated M1 | 0.12 / 1.46 | | 0.754 / 0.068 |
| no-split-sibling arms (7p, 7q, 14q, 22q; least arm leakage) | 2.61 | 0.41–16.6 | 0.151 |
| call set (1,187,033 TEST pairs called) | 1.15, Fisher p 3 × 10⁻¹¹ | | |
| regulator classes: RBP (17 regs, 30 hits) / other (1,596, 101) / TF (137, 0 hits) / kinase (51, 0) / ribosomal (16, 0 hits from 631 top-K pairs) | 3.56 / 3.21 / 6.67* / 42.9* / 5.43* | | *Haldane artefacts of empty cells |

### 5.6 HEK293T TEST (secondary line; 798 regulators, 8,373,505 pairs, K = 8,374)

| method | OR | perm p | note |
|---|---:|---:|---|
| selected | **0.49** (CI 0.02–14.1) | 0.878 | only 6 overlaps, concordance 0.33 |
| context layer | 5.06 (0.88–29.1) | 0.905 | |
| B4e / B5 / B2 / NC1 | 14.5 / 35.9 / 1.40 / 1.73 | 0.485 / 0.145 / 0.267 / 0.107 | |

## 6. TRAIN tuning path (HCT116 TRAIN; `results/tuning_log.tsv` compare rows, `results/tuning/*.json`)

Adoption rule (DESIGN §8.2, `benchmark.compare` `adopt_rule_8_2`): Δlog-OR > 0 in **both** TRAIN arm-folds,
paired two-way bootstrap 90 % lower bound > 0, and Δlog-CMH bootstrap 90 % upper bound > 0. Budget 30
configurations; 29 used (23 pre-registered incl. M1, 6 post-hoc D6). Greedy, one tier at a time.

| step | candidate (vs incumbent) | Δlog-OR | folds | boot lb90 | decision |
|---|---|---:|---|---:|---|
| baselines (M4) | M1 1.59 · B5 19.7 · B4e 7.09 · B3 6.99 · B2 5.77 · B4 5.05 · B1 1.26 · B6 0.97 · B0 0.98 · NC1 0.75 | (ORs) | | | — |
| tier 1 arm handling | focal | +0.68 | −0.09 / +1.52 | −0.41 | rejected (unstable) |
| | leave-one-arm-out | +0.04 | −0.00 / +0.10 | −0.06 | rejected |
| | LOCO CN-PCs (20) | −0.02 | −0.05 / +0.09 | −0.18 | rejected |
| tier 2 G1 (vs M1) | τ = 10 (w 0 or 0.5) | +0.25 | +0.21 / +0.36 | +0.12 | passes |
| | τ = 30 | +0.51 | +0.43 / +0.73 | +0.33 | passes |
| | **τ = 100, w = 0** | **+0.79** | +0.65 / +1.12 | +0.56 | **adopted** (w = 0.5 identical) → OR 3.50 |
| tier 3 G2 context expression (vs g1_t100_w0) | 4 options (threshold 1/2 × w 0/0.5) | +0.05 | e.g. +0.05 / +0.03 | ≈ −0.005 to −0.06 | rejected |
| tier 4 A3 (vs g1_t100_w0) | **ρ 0.95, s 2, w_co 0** | **+0.36** | +0.43 / +0.28 | +0.04 | **adopted** → OR **5.00** |
| | w_co 0.25 | +0.34 | +0.43 / +0.23 | +0.05 | passes, smaller gain |
| | w_co 0.5 | +0.28 | +0.33 / +0.21 | −0.01 | rejected |
| | ρ 0.98 | −0.05 | | −0.31 | rejected |
| | s 5 | 0.00 | | −0.48 | rejected |
| tier 5 target standardisation | | −0.03 | +0.09 / −0.18 | −0.26 | rejected |
| tier 6 hinge (loss slope) | | −0.11 | −0.25 / +0.10 | −0.60 | rejected |
| tier 7 transforms | CN floor 2⁻⁵ / RINT | −0.01 / −0.16 | | −0.07 / −0.50 | rejected |
| D6 post-hoc (label-free, added after seeing TRAIN) | I5 compositional TPM | +0.05 | | −0.07 | rejected |
| | I4 A3 ρ 0.90 | +0.21 | −0.04 / +0.51 | −0.36 | rejected (unstable) |
| | I3 shrunk-IV score | +0.18 (ΔCMH +0.42) | +0.10 / +0.29 | −0.07 | rejected |
| | I6 lineage-half replication | +0.04 | −0.01 / +0.14 | −0.12 | rejected |
| | I1 context fusion γ = 1 / **γ = 2** | +1.50 / **+1.57** (OR 24.1) | +1.37 / +1.82 (γ 2) | +1.30 | passes §8.2 but fails pre-declared C1 (exact C1 24.5, Δlog −0.02, CI −0.09 to 0.06) → **reported, not selected** (D7) |
| | I1 γ = 2 vs B5 | +0.20 | +0.04 / +0.41 | −0.21 | not better than hub prior |
| | I2 regulator-context weight | — | | | not run |

Selected method on HCT116 TRAIN (`results/eval/g1_t100_w0__a3_r095_s2_w0/hct116_TRAIN.json`; **200 permutations**):
OR 5.00 (jk CI 3.62–6.90), perm p 0.005 (= floor 1/201), degree-perm p 0.055, CMH 1.09 (0.98–1.21), AUROC 0.505
(333 regulators), sign concordance 0.561 vs c₀ 0.483, within-regulator null mean 0.494, sign p 0.005,
sign-concordant OR 6.19, 401 overlaps. HEK293T TRAIN OR 2.47. C2 (fusion on NC1) TRAIN 13.9.

TRAIN sensitivity runs vs selected (Δlog-OR): leave-HCT116-out −0.005, CN-outlier +0.005, TP53-chr17 +0.003,
S1 −0.10, hidden factors +0.13 (folds +0.23/−0.11), minus-G1 **−0.91**, focal +0.04 (folds −0.53/+0.50).

TRAIN diagnostics (reviews/ideas/M5a_ideas.md): within-regulator AUROC 0.507; ~94 % of log-enrichment is
regulator allocation; a G1 regulator's t row correlates with its arm-only (B6) row at median r = 0.87.

## 7. Runtime and memory (Apple M4 laptop, 10 cores, 32 GB, no GPU; Accelerate BLAS)

| fact | value | source |
|---|---|---|
| selected method, one full run (12,141 × 18,276 scored + 1 in-process NC1 pass for calls) | 67.2 s wall (timer starts **before** data loading, stops before writing scores.npz); engine time of the real pass 38.3 s | meta.json `wall_time_s`, `diagnostics.wall_engine_s` |
| context layer run | 82.7 s; engine 54.1 s | `results/runtime_table.tsv` |
| M1 | 51.6 s wall, 82 s end-to-end (process incl. writing output), `/usr/bin/time -l` footprint 6.4 GB | docs/MODEL.md §9 |
| model variants, wall s | 27.5 (NC1) / 38.7 (focal) … 151.4 (hinge); 51.5–90.7 for the runs in the final evaluation | `results/runtime_table.tsv` |
| baselines | 2.4–12.2 s | runtime_table |
| macOS ru_maxrss | 14.2 GB (selected), 20.2 GB (context layer), 13.3 GB (M1) | runtime_table `peak_rss_gb` |
| footprint (`/usr/bin/time -l`) for variants on the selected method | 7.3–8.4 GB | docs/MODEL.md §10 |
| score file size | 2.9 GB (`scores.npz`, uncompressed float32 t, score, alpha + bool call) | file size |
| data loading (`io.load_real`, warm cache) | ≈ 1 s (0.97 s with checksums) | recomputed 2026-10-28 |
| one benchmark evaluation (1,000 perms + 1,000 bootstraps) | 1.4–3.2 s | summary.md "runtime" `eval_wall_s` |
| full locked TEST evaluation (55 evaluations, 23 runs/controls, 2 lines) | 4.5 min | `results/final/LOCK` timestamps |
| test suite | 187 tests collected | `pytest --collect-only` (2026-10-28) |

## 8. Discrepancies (REPORT.md / RESEARCH_TREE.md / README.md vs code and result files)

1. **REPORT §1 runtime, "52–91 s across all model variants":** `results/runtime_table.tsv` model runs range from 27.5 s (NC1) / 38.7 s (focal) to 151.4 s (hinge); w_co variants 107 s, I6 103 s. 52–91 s is the range of the model runs *in the final evaluation* (focal+G1+A3 51.5 → hidden factors 90.7).
2. **REPORT §1 runtime, "plus about 30 s of data loading end to end":** `scripts/run_model.py` sets `t0` before `dio.load_real`, so the 67 s already includes loading; `load_real` takes ≈ 1 s (measured, warm cache). The ≈ 30 s gap between wall time and end-to-end time (M1: 51.6 vs 82 s, MODEL.md §9) lies after the timer stops — most plausibly writing the 2.9 GB `scores.npz` plus interpreter start-up (inferred, not measured).
3. **REPORT §1 finding 3, G1 TEST gain "+0.97 … consistent with TRAIN (+0.79)":** +0.79 is G1 added to M1 without A3 (tier 2). The TRAIN analogue of the TEST contrast (selected vs selected-minus-G1) is +0.91 (`tuning_log.tsv` compare `a3_r095_s2_w0` vs selected, −0.91). Both support G1; they are different contrasts.
4. **RESEARCH_TREE M4 "[!] M1 has direction information: concordance 0.581 vs independence c0 0.476":** `results/eval/M1/hct116_TRAIN.json` `sign.sign_perm_p` = 0.657 (within-regulator null); for `g1_t100_w0` it is 0.159. Only the A3 incumbent reaches p = 0.005. M1's "direction information" holds only against the pooled c₀, which REPORT itself treats as the weaker null.
5. **TRAIN p-values (tree M5/M5a "perm p 0.005 / ≤ 0.005", sign p 0.005):** come from 200-permutation evaluations (`results/eval/*/hct116_TRAIN.json` `n_perm` 200; `tuning_log.tsv` rows at 20:07–20:46), so 0.005 = 1/201 is the minimum attainable. The committed JSONs for M1, g1_t100_w0 and the selected run are the 200-permutation versions. REPORT §2.3's "1,000 permutations" applies to TEST.
6. **RESEARCH_TREE status, "scores ≈ 2.2×10⁸ trans pairs":** meta.json `n_pairs` 221,888,916 (all pairs) vs `n_trans_pairs` 209,943,680 (≈ 2.1×10⁸). REPORT words it correctly ("≈ 2.2×10⁸ pairs").
7. **REPORT §6 (CMH limitation bullet):** duplicated, garbled text "(B4e alone scores 3.8), (B4e alone scores 3.8 on TRAIN and 4.1 on TEST)". The values are right (TRAIN 3.77, TEST 4.12).
8. **REPORT §4 regulator classes** omit the `ribosomal` class in `results/final/summary.md` (16 regulators, 631 top-K pairs, 0 hits). This is an omission, not a contradiction.
9. **Terminology, DESIGN §6 and RESEARCH_TREE M1 / tree:** "83 % of HEK293T trans hits come from non-expressed *targets*" and "10 % of targets give > 90 % of hits" use PIE's meaning (target = **perturbed** gene). In the project's own convention, target = B, so read these as non-expressed *perturbations/regulators* (DATA_PROFILE §10 table). REPORT §1.6 phrases it correctly.
10. **docs/MODEL.md §6, "CN-only regulators (426) have t_cis = 0":** true, but `scores.npz` `cis_t` has 428 zeros; 2 more regulators have constant expression in S2, so their cis t is NaN and is set to 0. Effect: none (all fail G1).

Clarifications (consistent, but easy to misstate):
- A3 competes each regulator against **every** in-scope CN gene in its segment, including genes never perturbed and genes not in the regulator list. In the worked example, MDM4 takes ~14 % of the posterior and was never perturbed in HCT116. Segmentation uses LOCO-residualised CN, not raw CN (MODEL.md §6).
- The selected method's TEST C1 ran in the harness's "score (approximation)" mode (`summary.md` C1 `c1_mode`). For a CN-only score with no separate target factor, this is the intended shuffle.
- The 67 s run time includes one NC1 permutation pass, which is used only for the call set, not for the ranking.

## Correction D9 (2026-10-28) — supersedes fact 4 / any "direction is modestly right" wording
The within-regulator sign-concordance permutation test is **anti-conservative**: the negative control C2_g2
(context fusion on permuted CN) passes it on HCT116 TEST — concordance 0.585 vs within-regulator null mean 0.532,
pooled c₀ 0.514, p = 0.001, n = 354 (results/final/C2_g2/hct116_TEST.json → metrics.sign). NC1 does not
(0.364 vs 0.342, p 0.52, n 44). Therefore the selected method's 0.718 vs 0.647 (p 0.012, n 131) is **not
established evidence** of correct effect direction. Writers must say: "direction: not established (test failed
its negative control)". See RESEARCH_TREE.md D9, REPORT.md §1 finding 2.
