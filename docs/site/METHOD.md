# Copy-number trans effects in cancer cell lines: the method

> AI-generated project, not yet verified by a human. See the Disclaimer in the [README](../../README.md).
> Interactive version: [index.html](index.html) · Journey: [JOURNEY.md](JOURNEY.md) · Insights: [INSIGHTS.md](INSIGHTS.md)

**In one sentence:** DepMapTR asks which genes change expression when another gene's copy number changes; it takes
matched copy-number and RNA-seq profiles of 1,436 DepMap cancer cell lines, scores all 210 million pairs of genes on
different chromosomes for how well the copy number of a regulator gene A predicts the mRNA of a target gene B once
confounders are removed, and outputs a ranked list of regulator→target pairs with a predicted direction of change.

**The method in three sentences**

1. Natural copy-number variation across cell lines is treated as a weak, confounded experiment: B's mRNA is regressed
   on A's copy number, adjusting for B's own copy number, cell-line covariates (lineage, ploidy, culture medium and
   others) and the copy number of every chromosome arm *outside* A's chromosome.
2. Only regulators whose own mRNA clearly rises with their copy number are kept (the [cis gate](#step-s2)), and each
   signal is shared among genes that are always gained or lost together with A, favouring the one whose mRNA tracks
   its dosage best ([attribution](#step-s3)).
3. Each pair's score is the adjusted t-statistic × the gate (1 or 0) × the attribution share; pairs are ranked by its
   size, and the opposite of its sign is the predicted change in B when A is knocked down ([Step 4](#step-s4)).

**What it showed.** On held-out regulators the top 0.1% of ranked pairs were 3.85 times richer in CRISPRi Perturb-seq
hits than the rest, but this is not significant under the pre-registered test (p = 0.068) and not significantly better
than the same method run on shuffled copy number. The one replicated result, that the cis gate helps (+0.97 log-OR,
95% CI 0.30–1.64), suggests copy number says something about *which regulators* have detectable downstream effects; no
information about *which targets* respond was detected, and the predicted direction is not established
([Results](#results)).

| | |
|---|---|
| **Input** | DepMap copy number and RNA-seq (release 26Q3, inferred), 1,436 cell lines, plus lineage and culture metadata |
| **Output** | A score for each of 12,141 regulators × 18,276 targets on other chromosomes, and a predicted up/down change on knockdown |
| **Benchmark** | CRISPRi Perturb-seq in HCT116 cells (Arc Institute X-Atlas/Orion), regulators on 14 held-out chromosome arms |
| **Cost** | 67 s for all pairs on a 10-core laptop CPU, no GPU; 6.4–8.4 GB memory footprint (macOS max RSS 14 GB) |

## Terms used here

| term | meaning |
|---|---|
| copy number (CN) | DepMap's relative copy number of a gene: 1 = the line's average (its ploidy), 2 = twice as many copies. The model uses log2 of it, clipped to ¼×–8×. |
| cis / trans | *Cis*: a gene's own copy number acting on its own mRNA. *Trans*: gene A's copy number acting on gene B on a *different chromosome*. Only trans pairs are scored. |
| regulator A, target B | A is the gene whose dosage varies (and, in the benchmark, is knocked down); B is the gene whose mRNA responds. |
| chromosome arm | Each chromosome has a short (p) and long (q) arm. Cancers often gain or lose whole arms, so thousands of genes change copy number together. |
| CRISPRi Perturb-seq | A screen that knocks down one gene per cell (CRISPR interference) and reads every cell's transcriptome. A *hit* is a (knocked-down gene, measured gene) pair with adjusted p < 0.05. |
| odds ratio at top 0.1% (OR) | The main score: how much more often Perturb-seq hits occur among the 0.1% best-ranked pairs than among all other pairs. OR 1 = no better than chance. |
| t | A t-statistic for the A→B association after adjustment. \|t\| around 2 is borderline, above 5 strong. |
| Δlog-OR | Difference in log odds ratio between two rankings on the same pairs; +0.69 means the OR doubled. CI is a 95% interval from leaving out one held-out arm at a time. |
| per-regulator AUROC, CMH | Two "within one regulator" checks: AUROC is the chance that a hit target outranks a non-hit target of the same regulator (0.5 = chance); CMH is an odds ratio pooled within regulator × target-expression strata (1 = no within-regulator enrichment). |

<a id="differs"></a>
## How it differs from standard approaches

**The usual ways of finding trans effects either mistake co-gained chromosome arms for regulation, or remove genuine
dosage signal together with the noise.** This method keeps copy number as the upstream "cause" and adjusts only for
confounders that are not downstream of A: cell-line covariates and the measured copy number of other arms.

| approach | what goes wrong | what this method does instead |
|---|---|---|
| **Naive CN–expression correlation** (A's copy number vs B's mRNA) | Arms are gained together across chromosomes, and lineage drives both copy number and expression. GREB1 looks like a strong target of PPP1R15B (raw t = 6.1) but falls to t = 1.5 after adjustment. Held-out OR 1.14. | Adjusts for lineage, covariates and every arm outside A's chromosome ([Step 1](#step-s1)). Gain over naive on held-out data +1.22 log-OR, but the 95% interval (−0.48 to 2.92) includes zero. |
| **Co-expression** (A's mRNA vs B's mRNA) | No direction: shared programmes and B's own regulators all create correlation; it does not answer what a copy-number change does. Held-out OR 6.03, numerically higher than this method (difference not significant). | Anchors every association on A's copy number, upstream of expression, and requires that it moves A's own mRNA ([Step 2](#step-s2)). |
| **Trans-eQTL with expression-PC or PEER correction** | In tumours copy number drives the expression of 99% of well-expressed genes [F1], so top expression axes partly *are* copy-number effects; removing them removes genuine trans hubs [K1] and can create false positives [D2]. | Never uses expression factors; adjusts for *measured* confounders. **Not tested here:** expression-PC correction was never run, so its harm is taken from the literature. |

**Why it matters.** A copy-number association is useful as a prediction of "what happens if this gene's dosage
changes" only if it survives adjustment for what co-varies with dosage. The cis gate is the "strong instrument"
condition of Mendelian randomisation [G2, B2], and attribution is single-effect fine-mapping [W4].

## The pipeline

![The final method](../briefing/figures/pipeline.png)

Five steps turn two data matrices into a ranked list; a separate evaluation decides whether the list means anything.
Each step below says what it does, why it is needed, what goes wrong without it, and which
[turning point](JOURNEY.md#turning-points) produced it.

<a id="step-s1"></a>
### Step 1 · Measure the A→B association with confounders removed

For each regulator A and each target B on a different chromosome, one linear regression across the 1,436 lines: B's
mRNA on A's copy number, adjusting for B's own copy number (plus a flag for complete loss of B where at least 3 lines
have one), 56 covariate columns (lineage in 32 groups, ploidy, whole-genome doubling, MSI, TP53, sex, culture medium,
serum-free medium, growth pattern, RNA library type), a measure of overall aneuploidy, and the arm-level copy number of
every arm *not* on A's chromosome (leave-own-chromosome-out, LOCO).

    y_B = β·x_A + γ·[x_B, deletion_B] + Σ_(arms k not on chr(A)) δ_k·arm_k + covariates + error

- **Why:** cancers gain and lose arms in recurring combinations (88 co-occurring arm pairs [T1, S4]); without
  adjustment A's copy number stands in for those arms. A's own chromosome is excluded because adjusting for A's own arm
  would remove most of A's own variation.
- **Without it:** GREB1 appears to respond to PPP1R15B with raw t = 6.1; lineage and covariates bring it to 2.7, other
  arms' copy number to 1.5.
- **Evidence:** three alternatives (also adjusting A's own arm, leaving out only A's arm, copy-number PCs) did not beat
  it in training. On its own Step 1 is weak: training OR 1.59, permutation p 0.25.
- **Speed:** closed form (Frisch–Waugh–Lovell: residualise on the covariates once per chromosome, then partial
  correlations by matrix products); matches a standard regression package to 1e-8.
- Produced at [TP1](JOURNEY.md#tp1); confirmed at [TP6](JOURNEY.md#tp6).

<a id="step-s2"></a>
### Step 2 · Keep only regulators whose own mRNA follows their copy number (cis gate)

A passes if its own mRNA rises with its own copy number with adjusted t ≥ 10 (first-stage F ≥ 100, the usual measure of
instrument strength). Regulators that fail score 0 for every target. 5,403 of 12,141 pass.

- **Why:** copy number can only perturb A's targets through A if it changes A's product; weak instruments give biased,
  noisy answers [G2, B2].
- **Without it:** held-out OR falls from 3.85 to 1.46. The paired difference, +0.97 log-OR (95% CI 0.30–1.64), is the
  only gain that replicated from training (+0.79) to held-out data.
- **What it does not do:** it chooses regulators but does not reorder a regulator's targets (target-level statistic
  −0.05). It measures dosage response across all lineages, not in the screened line: MYCN passes through
  neuroblastoma amplicons but is not expressed in HCT116, and took 1,701 top-ranked training slots for 1 hit.
- Specified at [TP4](JOURNEY.md#tp4), adopted at [TP6](JOURNEY.md#tp6), replicated at [TP8](JOURNEY.md#tp8).

<a id="step-s3"></a>
### Step 3 · Share each signal among genes that are gained and lost together (attribution)

Neighbouring genes whose adjusted copy numbers correlate at r ≥ 0.95 with the running block average form blocks
(co-segments; 5,126 genome-wide). For each target B, assume one gene in the block drives B and compute each gene's
posterior share α: single-effect fine-mapping (SuSiE with one effect [W4]) with a closed-form approximate Bayes factor
(Wakefield's; no citation in the project's literature file).

    α(g | B) ∝ prior(g) × BF(g, B)     prior(g) ∝ exp(t_cis,g / 2)     prior effect sd σ₀ = 0.05

- **Why:** within an amplicon the trans t-statistics are nearly identical (7.0–7.5 for four genes in the example), so
  without a tie-breaker every passenger inherits the driver's score.
- **In practice:** the prior decides; α behaves almost like "give the block to the gene with the strongest cis
  response". Training gain +0.36 log-OR (lower bound +0.04), through choosing regulators, not targets.
- **Not validated:** not separately ablated on held-out arms.
- Specified at [TP4](JOURNEY.md#tp4), adopted at [TP6](JOURNEY.md#tp6), diagnosed at [TP7](JOURNEY.md#tp7).

<a id="step-s4"></a>
### Step 4 · Score, rank and predict a direction

S = t × gate weight (1 or 0) × α; pairs are ranked by |S|. The predicted change of B on knockdown of A is the opposite
of the sign of S (flipped again for the 3 regulators whose mRNA falls with more copies, cis t ≤ −3).

- **Why:** ranking needs one number; CRISPRi knockdown resembles copy loss.
- **Status: direction not established.** 71.8% of the 131 top-ranked held-out hits had the predicted direction,
  against 64.7% under a within-regulator shuffle (p = 0.012). A negative control built from shuffled copy number also
  passes this test (p = 0.001), so the test is miscalibrated.
- **Calls:** a 5%-FDR call set (|t| ≥ 3.016) has 10.7 million trans pairs and held-out OR 1.15; rank instead.
- Specified at [TP4](JOURNEY.md#tp4); status set at [TP8](JOURNEY.md#tp8).

<a id="step-s5"></a>
### Step 5 · Optional: re-rank for what one cell line can detect (context layer)

S_ctx = sign(t) · ρ_P · ρ_T², where ρ_P is the percentile of |S| and ρ_T is B's expression percentile in DepMap's own
HCT116 profile (no benchmark labels). Reported next to the method, not part of it.

- **Why:** Perturb-seq detects changes mainly in highly expressed genes (within a regulator, HCT116 target expression
  alone ranks hits with AUROC ≈ 0.70 in training).
- **What it adds:** held-out OR 18.4 (CI 10.9–30.9; p 0.015); beats the same fusion on shuffled copy number (+0.55,
  CI 0.07–1.03; post hoc, not multiplicity-corrected), which suggests the real copy-number evidence helps choose
  regulators; not different from a label-free hub baseline (22.8); shuffling targets within regulators barely
  changes it (17.3). It ranks by detectability, and was added after seeing training results.
- Produced at [TP7](JOURNEY.md#tp7).

<a id="step-e"></a>
### Evaluation · one sealed test on held-out chromosome arms

Universe: trans pairs whose regulator was knocked down in HCT116 and has variable copy number, and whose target is
measured in both datasets. Regulators split by whole arm into 28 training and 14 test arms (1,817 test regulators,
17.7 M test pairs); test labels sealed before tuning and opened once in a locked run.

- **Metric:** odds ratio of hits (adjusted p < 0.05) in the top 0.1% by |S|.
- **Significance:** regulator labels permuted 1,000 times, because hits cluster by regulator (the top 1% of
  perturbations give 55% of HCT116 trans hits); arm-jackknife confidence intervals.
- **Baselines and controls:** random, naive correlation, co-expression, target expression (variability; in HCT116), a
  hub prior (regulator cis strength × target expression), arm copy number alone; copy number shuffled within lineage;
  targets shuffled within each regulator; stratified odds ratio (CMH) and per-regulator AUROC.
- Shaped at [TP2](JOURNEY.md#tp2), [TP3](JOURNEY.md#tp3), [TP5](JOURNEY.md#tp5), [TP7](JOURNEY.md#tp7),
  [TP8](JOURNEY.md#tp8).

<a id="example"></a>
## Worked example: PPP1R15B on chromosome 1q

**One real regulator traced through every step with the project's own code** (numbers from
[`worked_example.json`](../explainer/worked_example.json)). PPP1R15B (CReP) directs protein phosphatase 1 to eIF2α; it
sits in a co-amplified block at 1q32.1 with MDM4. It is a *training* regulator, so only training labels are used.

1. **Copy-number variation.** Median relative CN 1.18 (1q is often gained); 45 lines with ≥ 2× gain, 1 with a 2×
   loss; focal amplicons up to 38.7× are clipped at 8×. HCT116 sits at 1.10.
2. **Cis gate.** Adjusted slope 1.04 log2 units per doubling, t = 23.9 (F = 572): passes. Block neighbours: PIK3C2B
   11.6 (passes), MDM4 20.4 (passes), LRRN2 8.5 (fails, scores 0).
3. **Trans association.**

   | target (chr) | raw r (naive t) | B's own CN only | covariates, no arms | **full model t** | effect β | training label |
   |---|---|---|---|---|---|---|
   | MYBL2 (20) | 0.20 (7.8) | 8.4 | 7.3 | **7.49** | +0.65 | hit, down |
   | NT5E (6) | −0.26 (−10.1) | −10.0 | −8.1 | **−7.35** | −1.18 | hit, up |
   | GREB1 (2) | 0.16 (6.1) | 5.9 | 2.7 | **1.46** | +0.13 | not a hit |

   Across all 16,388 trans targets, PPP1R15B's median |t| is 1.35, about twice noise (≈ 0.67): largely a 1q signature.
4. **Attribution** for target MYBL2:

   | gene | cis t | t with MYBL2 | α | knocked down in HCT116? | training hits |
   |---|---|---|---|---|---|
   | **PPP1R15B** | 23.9 | 7.49 | **0.934** | yes | 835 |
   | MDM4 | 20.4 | 7.34 | 0.066 | no | — |
   | PIK3C2B | 11.6 | 7.05 | 0.0002 | yes | 2 |
   | LRRN2 | 8.5 (fails) | 7.14 | 0.0001 | no | — |

   The data cannot separate the four genes; the cis prior decides. PIK3C2B (2 hits) gets almost no credit, PPP1R15B
   (835 hits) most of it: one anecdote, not proof. MDM4's share (6.6% here, 14% median over targets) cannot be checked.
5. **Score and direction.** MYBL2: S = 7.49 × 1 × 0.934 = 7.00; predicted down on knockdown; observed down. Top 10
   testable targets: POLD4, SLC39A13, **MYBL2** (hit, agrees), **NT5E** (hit, agrees), GNAI2, ZNF273, MVP, **RRAS**
   (hit, opposite), GSN, ZNF711.
6. **Benchmark.** 448 PPP1R15B pairs in the training top 0.1% (9th of 3,806 regulators). Hit rate: 11.8% for those
   pairs, 9.0% for its other targets, 0.22% for all training pairs. Choosing the regulator does most of the work;
   within-regulator odds ratio only 1.36. 34 of 53 top-ranked hits (64%) have the predicted direction; whether this
   beats chance is not established.
7. **Context layer.** The top 10 become highly expressed housekeeping-type genes (ACTG1, TPT1, MT2A, RPS12, PTMA,
   EEF1A1, TMSB10, ANXA2, RPS6, RPL7), with copy-number scores of only 2.5–6.6. Hits: 7 of 10 vs 3 of 10. What the
   layer adds is detectability, not biology.

## Pseudocode and where each step lives

n = 1,436 lines; A = 12,141 regulators; B = 18,276 targets. Paths under `src/depmaptr/`. Comments use the project's
internal ids: G1 = cis gate, A3 = attribution, NC1 = shuffled-copy-number permutation, I1 = context layer.

```text
 1  x[g] ← log2(clip(relCN[g], 1/4, 8))  for every CN gene g;   y[g] ← log2(TPM+1)
 2  Z ← lineage(32 lvls), ploidy, WGD, MSI, TP53, sex, medium, serum-free, growth, stranded (+ missing flags)
 3  arm[k] ← median x over genes on arm k  (41 arms; Xp/Xq sex-centred)
 4  for each regulator chromosome c:                                  # 23 groups (LOCO)
 5      D ← [1, Z, burden(arms ∉ c), arm[k] for k ∉ c]                # 98 columns
 6      Q ← pivoted-QR basis of D;  R(v) := v − Q Qᵀ v                # residualise once per group
 7      Ỹ ← R(y[B]) for all targets;  W̃_B ← R([x_B, 1(relCN_B<0.25) if ≥3 lines])
 8      for each CN gene g on c (all genes, needed by A3):
 9          t[g,B] ← partial t of R(x_g) vs Ỹ_B given W̃_B             # = OLS t of x_g in y_B~D+W_B+x_g
10          df[g,B] ← n − rank(D) − #active W columns − 1
11          t_cis[g] ← partial t of R(x_g) vs R(y_g)  (0 if g has no expression)
12      seg ← greedy_segments(R(x of all CN genes on c), ρ=0.95)     # join if r(gene, running mean) ≥ 0.95
13      for each segment s with m > 1 genes, each target B:
14          logpost[g] ← t_cis[g]/2 + ½log(1−λ) + ½λ·t[g,B]²,  λ = 0.05²/(0.05² + 1/df[g,B])
15          α[g,B] ← softmax_g(logpost)                               # SuSiE L=1; singleton α = 1
16  for each regulator A:
17      w_G1[A] ← 1 if t_cis[A] > 0 and t_cis[A]² ≥ 100 else 0        # G1: first-stage F ≥ 100
18      s_A ← −1 if t_cis[A] ≤ −3 else +1
19      for each target B on another chromosome:
20          S[A,B] ← t[A,B] · w_G1[A] · α[A,B]                         # final score; rank by |S|
21          kd_sign[A,B] ← −sign(S[A,B]) · s_A                         # predicted change of B on knockdown of A
22  # call set (reported, not used for ranking): one lineage-stratified permutation of R(x) rows (NC1);
23  t* ← min threshold with  #(|t_NC1| ≥ t*) / #(|t| ≥ t*) ≤ 0.05;   call[A,B] ← |t[A,B]| ≥ t*
24  # optional context layer (reported variant I1, γ=2; not the selected method):
25  ρ_P[A,B] ← percentile of |S[A,B]| among all finite trans entries (0 where S = 0)
26  ρ_T[B]   ← percentile of B's DepMap HCT116 log2(TPM+1) among all targets
27  S_ctx[A,B] ← sign(t[A,B]) · ρ_P[A,B] · ρ_T[B]²
28  # evaluation (benchmark.py): top 0.1% of trans universe pairs by |S| vs Perturb-seq p_adj < 0.05
```

| lines | step | file : function |
|---|---|---|
| 1 | copy-number transform | `io.py:cn_transform`, data via `io.py:load_real` |
| 2–3 | covariates; arm copy number | `covariates.py:build_Z`, `lineage_levels`, `arm_cn` |
| 5–6 | design and residualisation ([Step 1](#step-s1)) | `covariates.py:DesignBuilder.design`, `burden`; `assoc.py:Residualiser` |
| 7–10 | trans t and df ([Step 1](#step-s1)) | `pipeline.py:prepare`, `group_t` → `assoc.py:chunked_pair_partial`, `pair_partial` |
| 11, 17–18 | cis t, gate, sign rule ([Step 2](#step-s2)) | `assoc.py:diag_partial`; `gates.py:g1_pass`, `g1_weight`, `s_cis` |
| 12–15 | blocks and attribution ([Step 3](#step-s3)) | `attribution.py:greedy_segments`, `segment_alpha`, `log_abf` |
| 20–21 | score and direction ([Step 4](#step-s4)) | `pipeline.py:compose_score`; `benchmark.py:Evaluation.pred_kd_sign` |
| 22–23 | call set | `pipeline.py:nc1_permutation`; `calls.py:t_star`, `call_matrix` |
| 25–27 | context layer ([Step 5](#step-s5)) | `pipeline.py:context_fusion`, `trans_percentile`, `target_context_percentile` |
| 28 | evaluation ([E](#step-e)) | `benchmark.py:Evaluation.topk`, `compute_metrics`, `evaluate` |

<a id="results"></a>
## Results: what was shown, and how strongly

**On held-out regulators the copy-number method enriched Perturb-seq hits 3.85-fold, but this is not significant under
the pre-registered test, it does not beat simple expression-based baselines, and what enrichment there is appears to
come from choosing regulators rather than targets.** Held-out set: 14 arms in HCT116, 1,817 regulators, 17.7 M pairs, 34,343 hits.

| finding | evidence (held-out unless stated) | status |
|---|---|---|
| The cis gate improves enrichment. | OR 1.46 → 3.85; +0.97 log-OR (CI 0.30–1.64). Training: +0.79 when added, +0.91 when removed. | **Shown, replicated** |
| The ranking is enriched for hits. | OR 3.85 (CI 1.33–11.2); permutation p = 0.068; vs shuffled copy number +1.10 (CI −0.66 to 2.85); vs random +1.19 (CI 0.13–2.26). | **Suggestive, not significant** |
| It ranks targets within a regulator. | AUROC 0.500 (179 regulators); CMH 0.94 (CI 0.77–1.15); target shuffle 3.76 vs 3.85. | **No: none detected** |
| It predicts direction. | 71.8% vs 64.7% (p = 0.012, n = 131); shuffled-copy-number control also passes (p = 0.001); n.s. in 999-line set (p 0.13) and least-leaky arms (p 0.48). | **Not established** |
| It beats standard baselines. | Below co-expression 6.03, target expression in HCT116 7.39, hub prior 22.8 (no paired difference significant); above naive 1.14 (n.s.). | **Not shown** |
| It is robust to sample/covariate choices. | 999-line set 3.40; HCT116 removed 3.79; outliers removed 3.49; paralogs removed 3.82; hidden factors 2.90. | **Consistent** |
| It works in a second line. | HEK293T OR 0.49 (6 overlaps); not a DepMap line, 83% of hits from non-expressed perturbations. | **Not shown; benchmark weak** |

For scale: two CRISPRi screens agree at OR 181; a label-free hub prior reaches 22.8.

![Held-out TEST results](../briefing/figures/test_results.png)

![What signal there is concerns which regulators, not which targets](../briefing/figures/what_the_signal_is.png)

![Direction of predicted effects: not established](../briefing/figures/sign_concordance.png)

**Runtime.** One full run (2.2×10⁸ pairs plus one permutation pass) takes 67 s on a 10-core Apple M4 laptop CPU, 6.4–8.4
GB footprint (`/usr/bin/time`); baselines 2–12 s; one benchmark evaluation 1.4–3.2 s.

<a id="limits"></a>
## Limitations and when not to use it

**Assumptions.** Linear effects on the log scale, shared across lineages. Covariates and other arms' copy number capture
confounding; purity unavailable (lines assumed pure). A's copy number reaches B only through A or its block (attribution
handles neighbours only partly). DepMap release 26Q3 is inferred, not checksum-verified; ploidy/doubling from an older
release, imputed for 30% of lines.

**Not validated.** Target ranking (tested: none found); direction (test miscalibrated); any line other than HCT116;
attribution's held-out value; the harm of expression-factor correction (literature only). Nothing has been checked by a
human.

**When not to use it.**
- To pick targets for a regulator you care about: within a regulator the ranking is at chance.
- For regulators that fail the cis gate: their scores are 0 by construction.
- To name the causal gene inside a large amplicon: the prior decides, not the trans data.
- As a list of significant pairs: the 5%-FDR set has 10.7 M pairs (genomic inflation λ ≈ 2.1 from real arm structure).
- To predict the sign of a knockdown effect.

**Precision.** Only 14 held-out arms (CI 1.33–11.2). Secondary metrics are not multiplicity-corrected; only the
pre-registered primary metric is confirmatory. The context layer and its controls were added after seeing training
results.

<a id="credit"></a>
## Ideas borrowed from the literature

| idea | used for | source (keys from [LITERATURE.md](../../LITERATURE.md)) |
|---|---|---|
| Copy number drives most expression; expression factors absorb it | no expression PCs | [F1] Fehrmann 2015; [K1] Kang 2008; [D2] Dahl 2019 |
| Arm events are frequent and co-selected | other-arm adjustment; arm hold-out | [T1] Taylor 2018; [S4] Shukla 2020; [D3] Davoli 2013; LOCO from GWAS practice (no citation) |
| Closed-form residualisation at scale | engine | [S5] Matrix eQTL; [T2] tensorQTL |
| Copy number as an instrument; weak-instrument bias | cis gate | [G2] Goh 2013; [B2] Burgess & Thompson 2011 |
| Cis-anchored trans effects agree with Perturb-seq | motivation | [W1] Warmerdam 2026 (OR 12.4 on 41 pairs) |
| Single-effect fine-mapping | attribution | [W4] SuSiE; closed-form Wakefield Bayes factor from a reviewer agent (not in LITERATURE.md) |
| Co-expression does not track perturbation effects | co-expression demoted | [A1] Aguirre 2025 |
| Supervised hidden factors | sensitivity analysis | [L2] GBAT; [F2] PANAMA |
