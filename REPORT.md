# REPORT — DepMapTR: predicting trans-transcriptomic consequences of copy-number change

Date: 2026-10-27 · Repository: this directory (git) · Brief: PROMPT.md · Design: DESIGN.md (v3.1 + decision
notes D3–D8a in RESEARCH_TREE.md) · All numbers below come from committed files in `results/` (TEST:
`results/final/`, one LOCKed run; TRAIN: `results/tuning/`, `results/m4/`, `results/tuning_log.tsv`).

## 1. Summary

We built a closed-form, genome-wide tool that estimates, for every trans gene pair (A on a different chromosome
from B), the effect of A's copy number (CN) on B's mRNA across 1,436 DepMap cell lines. The model adjusts for
B's own (cis) CN, arm-level CN of all other chromosomes, genome-burden, lineage, ploidy/WGD, MSI, TP53, sex,
media and RNA batch. It then gates regulators on a strong cis CN→mRNA instrument and attributes signal among
co-segmented genes. We evaluated it against Arc Institute PIE X-Atlas/Orion CRISPRi Perturb-seq DE tables on
regulators held out by chromosome arm, which were sealed outside the repository before any tuning.

**Held-out (TEST) results.** HCT116 is the primary line; there are 14 held-out arms, 1,817 regulators and
17.7 M universe pairs.

| method | OR @ top 0.1% | jackknife 95% CI | regulator-label permutation p | sign-concordant OR | concordance (c₀; within-regulator null mean) | sign perm p | score time (s) |
|---|---:|---|---:|---:|---|---:|---:|
| **Selected CN-only tool** (M1 + G1 + A3) | **3.85** | 1.33–11.2 | 0.068 | 5.84 | **0.718** (0.541; 0.647), n=131 | **0.012** | 67.2 |
| Selected minus cis gate (M1 + A3) | 1.46 | 0.91–2.33 | 0.23 | 1.89 | 0.72 (0.49), n=50 | 0.006 | 66.1 |
| Context layer (selected + HCT116 target-detectability fusion) | **18.4** | 10.9–30.9 | **0.015** | 21.1 | 0.549 (0.507) | 0.22 | 82.7 |
| B5: DepMap-only hub prior (regulator cis×CN-SD × target HCT116 expression) | 22.8 | 3.3–159 | 0.043 | 26.2 | = c₀ (no sign information) | — | 2.4 |
| B4e: target expression in DepMap HCT116 | 7.39 | 5.55–9.83 | 0.53 | — | — | — | — |
| B4: target expression SD | 5.73 | 4.23–7.76 | 0.53 | — | — | — | — |
| B2: co-expression | 6.03 | 2.7–13.4 | 0.095 | 7.83 | 0.843 (0.845) | 1.0 | 4.8 |
| B0 random | 1.17 | 0.96–1.41 | 0.15 | 1.25 | — | — | 2.7 |
| NC1 permuted CN | 1.28 | 0.48–3.40 | 0.16 | 0.83 | — | — | 27.5 |

The findings:
1. **The CN-only tool's enrichment carries over to TEST in size and direction, but it is not significant under
   the pre-registered test.**
   - OR 3.85 on TEST. On TRAIN, M1 reached 1.59 and the selected tool 5.0. The TEST ablation without G1
     reaches 1.46.
   - The jackknife CI excludes 1, and the paired gain over random scores is significant (Δlog-OR +1.19, CI
     0.13–2.26).
   - The pre-registered regulator-label permutation gives p = 0.068, and the degree-stratified permutation
     p = 0.061. The paired gain over the permuted-CN control NC1 is not significant (+1.10, CI −0.66 to
     2.85). Fisher p is 2×10⁻³⁶, but only under an independence assumption we reject.
   - Target-detectability baselines (B4e 7.39, B4 5.73) and co-expression (B2 6.03) all score above the
     CN-only tool on pooled OR. None of these paired differences is significant.
2. **The predicted direction of effects looks better than chance, but the test is not trustworthy** (correction
   D9, 2026-10-28).
   - Among the 131 overlapping pairs, 71.8% are concordant (CN loss ↔ knockdown). The pooled independence
     expectation is 54.1%, and the within-regulator permutation null mean is 64.7% (p = 0.012).
   - **However, a negative control also passes the same test.** C2 builds the same scores from permuted CN:
     it gives 58.5% vs a null mean of 53.2%, p = 0.001, n = 354.
   - The within-regulator sign permutation ignores the fact that targets respond in correlated modules, so it
     is anti-conservative. Its calibration was never checked on negative controls (only the
     enrichment-permutation test was; DESIGN §6).
   - The sign result is therefore **not established**. It is also not robust in S1 (p = 0.13) or the
     no-split-sibling subset (p = 0.48), and the overlaps come from a few regulators (6 with ≥ 5 pairs).
3. **The cis-instrument gate (G1) is the key ingredient.** On TEST it improves log-OR by +0.97 (paired
   jackknife 95% CI 0.30–1.64), consistent with the matching TRAIN ablation (+0.91; the tier-2 adoption gain was +0.79).
4. **The tool has no detectable information for ranking targets within a regulator.**
   - Per-regulator AUROC is 0.500 (179 regulators), and the regulator-stratified CMH OR is 0.94 (CI
     0.77–1.15, which bounds any lift at about 15%).
   - The pair-shuffle control C1 gives 3.76 vs 3.85 (Δlog-OR +0.02). Its jackknife CI (−0.75 to 0.79) is
     too wide for C1 on its own to be informative.
   - DepMap CN tells us (weakly) *which* regulators have trans consequences detectable by CRISPRi. It does not
     tell us *which* of their targets respond, and the evidence on the *direction* is not established (finding 2).
5. **The context layer adds target detectability and reaches OR 18.4.**
   - It is not significantly different from the DepMap-only hub prior B5 (Δlog-OR −0.22, CI −2.29 to 1.86).
   - It does beat the same fusion built on permuted CN (C2; +0.55, CI 0.07–1.03). So the real CN evidence
     improves regulator allocation; this is real t plus A3, since C2 has neither.
   - After removing the top 1% hub regulators and targets, it reaches 14.8 versus B5's 11.8 (no paired
     test).
6. **HEK293T (secondary) does not replicate for the CN-only tool** (TEST OR 0.49; only 6 pairs overlap). The
   context layer reaches 5.1 there (perm p 0.91), below B4e (14.5) and B5 (35.9). HEK293T is not a DepMap
   line, and 83% of its trans hits come from non-expressed (likely artefactual) perturbations
   (docs/DATA_PROFILE.md §11).

**Runtime** (10-core laptop CPU, no GPU; per-method score times are in the table above):
- One full tool run scores 12,141 regulators × 18,276 targets ≈ 2.2×10⁸ pairs, including the NC1 permutation
  pass used for the call set. Compute time is 67 s for the selected method and 83 s for the context layer
  (27–151 s across all model and control runs), wall-clock including data loading (about 1 s).
- Peak memory footprint (`/usr/bin/time -l`) is about 6.4 GB for M1 and 7.3–8.4 GB for variants built on the
  selected method. macOS `ru_maxrss` reports 14 GB for the CN-only tool and 20 GB for the context layer.
- Evaluating one score file on the benchmark takes 1.4–3.2 s.
- The whole final TEST evaluation (55 evaluations of 23 runs over 2 lines, 1,000 permutations each) took
  4.5 min.

See `results/runtime_table.tsv` and the runtime section of `results/final/summary.md`.

## 2. Approach

### 2.1 Data (DATA.md has sources, versions, checksums and licences)
- **DepMap Public 26Q3.** The release is inferred, not md5-verified, because the portal file listing is behind
  a bot wall; DATA.md lists the evidence. We used `OmicsCNGeneWGS` (linear relative CN, ploidy-normalised)
  and `OmicsExpressionTPMLogp1HumanProteinCodingGenes`, with the default profile per model.
- **Sample set S2:** n = 1,436 cell lines. Organoids, engineered models and duplicate patients are excluded.
  The PureCN-complete set (S1, n = 999 within S2) was rejected as primary because it drops HCT116 itself and
  45% of MSI-H lines (data audit E2).
- **Covariates:**
  - PureCN ploidy/WGD and MSIScore are from 24Q4, the last release with them; missing values get
    lineage-median imputation plus an indicator.
  - TP53 status is from the 24Q4 damaging/hotspot mutation matrices.
  - RNA `Stranded`. Sequencing platform was dropped because it is 93% missing (D3).
  - Purity is unavailable; cell lines are assumed ≈ 1.
- **Benchmark:** PIE X-Atlas/Orion `de/*.parquet` (HF commit `abfa3da`), the pre-computed Mann–Whitney DE
  tables only. No raw Perturb-seq data were used.

### 2.2 Model (DESIGN §4–§5; code in `src/depmaptr/`)
- **Core (M1, LOCO reduced form):**
  y_B = β·x_A + γ·[x_B, 1(homdel_B)] + Σ_{arms k ∉ chr(A)} δ_k·arm_k + Z·η + ε.
  - x is log2 relative CN clipped to [2⁻², 2³].
  - Z holds the shared covariates: lineage (32 levels), ploidy, WGD, MSI, TP53, sex, base medium,
    serum-free, growth pattern, stranded, and genome burden computed leave-one-chromosome-out on autosomes.
  - It is fitted by Frisch–Waugh–Lovell, with one pivoted-QR residualisation per regulator chromosome and
    closed-form partial correlations to t. It matches statsmodels OLS to 1e-8 on 480 real pairs.
- **Why these covariates:**
  - lineage, media and batch are the strongest expression confounders;
  - arm CN of the *other* chromosomes absorbs co-selected arm events, the main trans confounder;
  - B's cis CN is required, and a homozygous-deletion indicator handles nonlinearity;
  - LOCO burden, ploidy and WGD capture the generic aneuploidy-stress axis;
  - MSI and TP53 are upstream drivers of CN burden and of distinct programmes.
  - Expression PCs/PEER were deliberately *not* used, because they absorb genuine CN-driven trans hubs. A
    regulator-conditional factor sensitivity analysis is reported in §4.
- **G1 cis-instrument gate:** a regulator is kept only if its own CN→mRNA first-stage F ≥ 100, i.e.
  t_cis ≥ 10. Others get weight 0.
- **A3 attribution:** a SuSiE L=1 single-effect posterior over the genes of each co-segment, using greedy
  segments at r ≥ 0.95. The prior is ∝ exp(t_cis/2); σ₀ = 0.05 is fixed a priori.
- **Score and sign:** S_AB = t_AB·w_G1·α_{A|B}. The predicted knockdown sign is −sign(S)·s_A, where s_A = −1
  only if t_cis ≤ −3.
- **Call set (MVP rule):** |t| ≥ t* such that the NC1-permutation FDR ≤ 5%.
- **Context layer (reported variant):** S = sign(t)·ρ_P(|S|)·ρ_T(B)². ρ_P is the percentile of the CN-only
  score and ρ_T is B's expression percentile in DepMap HCT116. It is label-free, but it ranks consequences by
  their *detectability* in the context line rather than by effect size.

### 2.3 Evaluation (DESIGN §6; `src/depmaptr/benchmark.py`, reviewed independently)
- **Universe (per line):**
  - A: perturbed in that line, mapped 1:1 to DepMap, autosomal or X, with variable CN (≥ 20 lines with
    |log2 relCN| ≥ 1; X is sex-centred);
  - B: consistently tested in Perturb-seq (control mean ≥ 0.05) and measured in DepMap;
  - chrom(A) ≠ chrom(B).
- **Positives:** PIE p_adj < 0.05 (BH within perturbation). Sensitivity sets are p_adj < 0.01 and
  |Ueffect| ≥ 0.2 (D4).
- **Primary metric:** Fisher OR among the top 0.1% of universe pairs ranked by |score|.
- **Significance:** regulator-label permutation (1,000 permutations), because pairs sharing a regulator are not
  independent. CIs come from a delete-one-arm jackknife.
- **Sign:** concordant iff the predicted knockdown sign equals sign(LFC). The test is against the independence
  expectation c₀, using within-regulator sign permutation.
- **Held-out split:** whole chromosome arms (triplets sorted by regulator count, seed 20261027), giving 14 TEST
  and 28 TRAIN arms.
  - TEST labels and the raw DE tables were physically moved out of the tree and made unreadable, with their
    hashes committed, before any evaluation.
  - One final run: its LOCK was committed before any TEST file was read (commits 603a7f5 → e9085c7).

## 3. How the method was chosen (TRAIN only; 29 of 30 budgeted configurations including M1; the 7 D8 sensitivity/ablation runs were evaluated on TRAIN for reporting only)

Greedy pre-registered tiers on HCT116 TRAIN were adopted only if the log-OR gain was positive in both TRAIN
arm-folds **and** the paired-bootstrap 90% lower bound was > 0 (§8.2). Full log: `results/tuning_log.tsv`.

| tier | options (Δlog-OR vs incumbent) | decision |
|---|---|---|
| baselines (M4) | M1 1.59; B5 19.7; B4e 7.1; B2 5.8; B0 0.98; NC1 0.75 | — |
| 1 arm handling | focal +0.68 (folds −0.09/+1.52), leave-one-arm-out +0.04, CN-PCs −0.02 | none adopted |
| 2 G1 cis gate | τ=10 +0.25, τ=30 +0.51, **τ=100 +0.79** | **adopted** |
| 3 G2 context expression of A | +0.05 | not adopted (redundant after G1) |
| 4 A3 attribution | **ρ 0.95/s 2/w_co 0 +0.36**; w_co 0.25 +0.34; ρ 0.98 −0.05; s 5 0.00 | **adopted** |
| 5–7 | target standardisation −0.03; hinge (loss slope) −0.11; CN floor 2⁻⁵ −0.01; RINT −0.16 | none |
| D6 post-hoc (label-free) | compositional TPM correction +0.05; ρ 0.90 +0.21 (unstable); shrunk-IV score +0.18 (ΔCMH +0.42); lineage-replicated t +0.04 | none |
| D6 context layer | detectability fusion γ=1 +1.50, **γ=2 +1.57** | passes §8.2 but fails pre-declared pair-shuffle control C1 → **reported, not selected (D7)** |

Key TRAIN diagnostics (reviews/ideas/M5a_ideas.md):
- Within a regulator, |t| ranks targets at noise level (AUROC 0.507).
- About 94% of log-enrichment comes from choosing regulators.
- Score rows are dominated by arm-level trans signatures (median r = 0.87 with the arm-only score).

## 4. Further TEST results

All of these are in `results/final/summary.md`.

**Robustness of the selected tool (HCT116 TEST OR@0.1%):**
| analysis | OR | note |
|---|---:|---|
| top 0.01% / 1% | 2.92 / 4.77 | |
| positives p_adj < 0.01 / \|Ueffect\| ≥ 0.2 | 4.05 / 2.65 | |
| paralog pairs excluded | 3.82 | |
| hub-excluded (top 1% regulators and targets removed) | 2.69 | Fisher p 7×10⁻⁶ |
| S1 sample set / leave-HCT116-out | 3.40 / 3.79 | |
| CN-outlier lines excluded / TP53 dropped for chr17 | 3.49 / 3.73 | |
| regulator-conditional hidden factors | 2.90 | |
| focal arm handling | 4.50 | Δ vs selected not significant |
| regulators whose own gene is knocked down (class i) | 1.61 | |
| regulators not expressed (class iii, empirical null) | 0.23 | 1 hit |
| NC2 (not expressed in benchmark or DepMap) | 0.12 (selected), 1.46 (ungated M1) | few hits; consistent with null |
| no-split-sibling regulators (other arm also in TEST) | 2.61 | |
| call set (M1 \|t\| ≥ t*, NC1-FDR 5%; 1.19 M TEST-universe pairs, 10.7 M genome-wide) | 1.15 | Fisher p 3×10⁻¹¹; identical for selected, minus-G1 and context layer, because t* thresholds the M1 t |

- **Leaving HCT116 itself out of DepMap changes nothing** (Δlog-OR −0.02 vs selected, CI −0.08 to 0.05), so the result is not driven by the
  benchmark line being in the training panel.
- **PC2, the empirical ceiling:** HEK293T Perturb-seq hits used to predict HCT116 hits for shared regulators
  give OR 181. This is the scale of agreement between two CRISPRi screens, far above anything observational.
- **Regulator classes** (HGNC groups, rules fixed before results):
  - RBP OR 3.56 (30 hits; Fisher p 1×10⁻⁸);
  - ribosomal: 0 hits among 631 top-K pairs (16 regulators);
  - "other" 3.21 (101 hits);
  - TF and kinase have **0 hits** in the top-K. Their displayed ORs (6.7, 42.9) are Haldane-corrected
    artefacts of empty cells and are uninformative.
  - The CN signal is not TF-enriched (0 TF hits, although power is low), consistent with the brief.

## 5. Interpretation

DepMap copy-number variation acts as a weak and heavily confounded natural perturbation; the evidence that it is informative at all is suggestive rather than conclusive on TEST (§1):
- **Which regulators.** Regulators whose own mRNA strongly tracks their CN are the ones whose CRISPRi
  knockdown has detectable trans consequences. Restricting to them roughly doubles to triples enrichment.
- **Direction (not established).** 72% of overlapping pairs are concordant, against a 65% within-regulator
  null. But a permuted-CN negative control passes the same test, so the test is miscalibrated and this is not
  evidence (D9).
- **Not which targets.** Within a regulator, which targets respond in a single-cell CRISPRi screen of one line
  is not detectably predicted by cross-line CN associations. Perturb-seq power is dominated by target expression
  depth, and DepMap trans associations are dominated by arm-level co-dosage signatures.
- **The pooled benchmark rewards hubness and detectability.** A label-free DepMap prior (B5) reaches OR ≈ 20.
  An honest evaluation of observational trans inference must therefore report hub baselines and
  within-regulator controls next to the headline OR, as we do.

## 6. Limitations

- **TEST precision is limited by 14 held-out arms.** The CIs are wide (selected 1.3–11.2), and the
  regulator-label permutation p is 0.068 for the CN-only tool.
- **Several choices were post hoc.** D6 added label-free variants and the C1/C2 controls after seeing TRAIN
  results; they are documented and were capped by the 30-configuration budget. The adoption rule is strict,
  so gains with one unstable fold were rejected; some may be real (focal, I4).
- **HEK293T benchmark quality is low.** It is not a DepMap line, and its hits are dominated by non-expressed
  perturbations. The CN-only tool fails there on TEST, with very few overlaps.
- **The release is unverified** (26Q3 by inference). PureCN covariates are frozen at 24Q4 and missing for 30% of S2 lines (437/1,436; imputed, with an indicator). Purity is unavailable.
- **The LOCO model attributes each arm-level trans signature to every gene on the arm.** A3 attribution
  helps regulator allocation, but cannot identify the causal gene without extra information.
- **The CMH stratification is coarse.** The frozen CMH (pan-DepMap expression deciles) does not remove
  HCT116 detectability (B4e alone scores 3.8 on TRAIN and 4.1 on TEST), so CMH was supplemented by the C1 control and the D6
  HCT116-stratified CMH.
- **The call set is large** (10.7 M trans pairs genome-wide; 1.2 M in the TEST universe) because λ_GC ≈ 2.1 reflects genuine arm-level structure. Its
  enrichment is modest (OR 1.15). For use, ranking is preferable to calling.

- **The sign-concordance test is miscalibrated.** The permuted-CN control C2 passes it (p = 0.001), so the
  selected method's p = 0.012 is not evidence of correct directions. It is also not robust (S1 p = 0.13;
  no-split-sibling p = 0.48). A calibrated test would need to preserve target-module structure, for example
  by permuting regulators rather than target signs.
- **The CN-only tool does not significantly beat NC1 in the paired test**, and its pooled OR is below B2, B4
  and B4e.
- **The no-split-sibling subset is weaker** (OR 2.61, CI 0.41–16.6, perm p 0.15). It contains regulators
  whose chromosome's other arm is also in TEST, so it has the least arm-level leakage from TRAIN.
- **Many secondary TEST metrics are reported without multiplicity adjustment.** Only the pre-registered
  primary metric should be read confirmatorily.

## 7. Follow-ups and remaining ideas

1. **Allele-specific / focal-event instruments.** Use focal amplicons and homozygous deletions, rather than arm
   dosage, to break arm co-dosage, for example by estimating effects only from lines with focal events at A.
2. **Proteomics mediation.** DepMap/CCLE proteomics would test whether A's protein, not only its mRNA, tracks
   CN (dosage compensation); a protein-level instrument gate could sharpen G1.
3. **Context-specific estimation.** Weight lines by similarity to the benchmark line, or test
   lineage-interaction models. Pre-specified heterogeneity-aware pooling was not run.
4. **Multi-regulator (joint) models.** Use sparse multivariate regression over segment-level CN, so that
   co-gained regulators compete per target.
5. **Better benchmarks for pair-specific evaluation.** Within-regulator metrics with fine detectability
   strata, and benchmarks restricted to regulators with confirmed knockdown and expressed targets.
6. **Revisit near-misses with a fresh hold-out.** These include focal arm handling (TEST 4.50 vs 3.85), the
   shrunk-IV score (TRAIN ΔCMH +0.42) and coarser segments.
7. **Verify the DepMap release md5s manually** in a browser (DATA.md).

## 8. Reproducing

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.lock.txt && .venv/bin/pip install -e .
scripts/data/run_all.sh        # data foundation: needs starting_data/ + PIE DE parquet (HF abfa3da) in data/external/pie_x_atlas_orion/
.venv/bin/python scripts/data/download_paralogs.py && .venv/bin/python scripts/data/build_regulator_classes.py
.venv/bin/python scripts/benchmark/build_benchmark.py     # benchmark tables + split from the PIE DE tables; fresh clone only (it regenerates TEST labels, which this repo keeps sealed)
.venv/bin/python scripts/run_model.py --variant g1_t100_w0+a3_r095_s2_w0 --out data/results/scores/g1_t100_w0__a3_r095_s2_w0
.venv/bin/python scripts/evaluate.py g1_t100_w0__a3_r095_s2_w0 --line both   # TRAIN only
.venv/bin/python -m pytest -q tests                       # 187 tests (SIM known-answer, OLS agreement, metrics)
.venv/bin/python scripts/runtime_table.py                 # results/runtime_table.tsv
# other score runs: every run id in results/final/selected.yaml, via run_model.py --variant <spec with '+'>
```

Seeds are in `config/seeds.yaml` (mirrors `results/variants.yaml`: [20261027, 1, 2, 3, 4]; streams use `SeedSequence(20261027)` spawns logged per run). The TEST run
is locked (`results/final/LOCK`), and `scripts/final_eval.py` refuses to run again while the LOCK exists.
