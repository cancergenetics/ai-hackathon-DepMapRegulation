# Worked example: one gene through the DepMapTR method

This page follows one real gene, **PPP1R15B**, through every step of the final ("selected") DepMapTR method and
then compares its predictions with the CRISPRi Perturb-seq benchmark. Every number comes from
`worked_example.json`. That file is written by `make_worked_example.py` in this folder, which re-runs the
project's own functions on the real data:

```bash
.venv/bin/python docs/explainer/make_worked_example.py      # about 20 s on a laptop
```

The script also checks its recomputed values against the committed score file of the selected run. The
largest differences are 4 × 10⁻⁷ in t and 3 × 10⁻⁸ in α, which is float32 rounding.

PPP1R15B is a **TRAIN** regulator (chromosome arm 1q is a tuning arm). Only TRAIN Perturb-seq labels are used
here; the held-out TEST labels stay sealed.

## The question the tool answers

The question: *if a cell line carries more or fewer copies of gene A, which other genes (B, on other
chromosomes) go up or down in expression?*

The tool learns this from natural copy-number (CN) variation across 1,436 DepMap cancer cell lines. It is
judged against a CRISPRi Perturb-seq screen in HCT116 cells. In that screen A is knocked down (an artificial
"copy loss"), and the screen records which B genes change.

**Why PPP1R15B?**
- It has a strong cis link: its own mRNA tracks its copy number very closely.
- It sits in a multi-gene co-amplified block on 1q32 that also contains the oncogene **MDM4**, so the
  attribution step has real work to do.
- Its knockdown in HCT116 has many trans effects: 835 significant TRAIN hits, which ranks it 28th of 3,806
  TRAIN regulators.

Biologically, PPP1R15B (also called CReP) targets protein phosphatase 1 to eIF2α. Losing it switches on the
integrated stress response, so broad downstream effects are expected.

## Terms used

| term | meaning |
|---|---|
| relCN | DepMap relative copy number. 1 means the same copy number as the line's average (its ploidy); 2 means twice as many copies |
| x | the CN value the model uses: x = log2(relCN), clipped to [−2, 3] (¼× to 8×). x = 1 means 2× gain; x = −1 means 2× loss |
| y | expression, log2(TPM + 1) |
| covariates | things that differ between cell lines and could create false associations: tissue lineage (32 groups), ploidy, whole-genome doubling, MSI status, TP53 status, sex, culture medium, growth pattern, RNA library type, overall aneuploidy ("burden"), and **arm-level CN of every chromosome arm except those on A's own chromosome** |
| t | a t-statistic: how strongly two things are associated after removing the covariates. \|t\| ≈ 2 is borderline; \|t\| > 5 is strong |
| cis | A's own CN → A's own mRNA |
| trans | A's CN → mRNA of a gene B on a different chromosome |
| TRAIN hit | a pair (A, B) where knocking down A changed B in HCT116 Perturb-seq (adjusted p < 0.05), in the tuning half of the benchmark |

## Step 1 — How much does PPP1R15B's copy number vary across 1,436 cell lines?

PPP1R15B sits at 204.4 Mb on 1q (1q32.1). Arm 1q is gained in many cancers, and the distribution shows it.

| quantile | min | 1 % | 5 % | 25 % | median | 75 % | 95 % | 99 % | max |
|---|---|---|---|---|---|---|---|---|---|
| relCN | 0.40 | 0.74 | 0.85 | 1.05 | **1.18** | 1.41 | 1.81 | 2.49 | 38.7 |
| x = log2 | −1.31 | −0.44 | −0.23 | 0.07 | 0.23 | 0.49 | 0.86 | 1.31 | 3.0 (clipped) |

- The typical line carries about 18 % more copies than its ploidy.
- 45 lines have at least a 2× gain (x ≥ 1); only 1 has a 2× loss.
- That gives 46 lines with |x| ≥ 1, above the 20 needed to count as a "variable-CN" regulator.
- **Highest:** CCLF_NEURO_0014_T (brain, relCN 38.7), CCLF_NEURO_0007_T (brain, 14.5) and RMZ (soft tissue,
  11.9). These are focal high-level amplifications, consistent with the known 1q32/MDM4 amplicon. The clip at
  8× (x = 3) stops these few lines from dominating the fit.
- **Lowest:** JHOS-2 (ovary, 0.40).
- **HCT116 itself:** relCN 1.10 (x = 0.14), close to neutral.

## Step 2 — Is this a usable instrument? (own expression vs own CN, gate G1)

Before trusting any trans effect, the tool asks whether A's copy number actually changes A's own mRNA. If it
does not, "more copies of A" is not a perturbation of A at all.

The tool fits y_A against x_A while holding the covariates fixed: 98 design columns, 1,337 degrees of freedom.

| | raw (no covariates) | adjusted (the model) |
|---|---|---|
| correlation | r = 0.55 | partial r = 0.55 |
| slope (log2 mRNA per log2 CN) | 1.08 | **1.04** |
| t | — | **t_cis = 23.9** (first-stage F = t² = 572) |

- A slope of about 1 means doubling the gene roughly doubles its mRNA, i.e. proportional dosage.
- A statsmodels refit of the full regression gives the same t (23.92) and slope (1.041).
- **Gate G1:** a regulator passes only if t_cis > 0 and t_cis² ≥ 100, i.e. t_cis ≥ 10. PPP1R15B passes easily,
  so its weight is w_G1 = 1. A regulator that fails gets weight 0: its whole row of scores becomes 0.
  Genome-wide, 5,403 of 12,141 regulators pass.
- **Sign rule:** s_A = −1 only if t_cis ≤ −3 (a gene whose mRNA goes *down* with more copies). Otherwise
  s_A = +1. Here s_A = +1.

`scatter_subsample` in the JSON holds about 400 (x_A, y_A) points for plotting; HCT116 is the last point.

## Step 3 — Trans association: PPP1R15B's CN vs other genes' expression

For every gene B on another chromosome, the tool fits:

> y_B ~ x_A + (B's own CN, plus a homozygous-deletion flag for B if ≥ 3 lines have one) + all covariates
> including arm CN of every arm not on chromosome 1

The coefficient on x_A is the estimated effect. Its t-statistic is the evidence. The engine does this for all
18,276 targets at once with a closed-form shortcut (Frisch–Waugh–Lovell), not one regression per pair. The
table refits three pairs one at a time with statsmodels to show what the adjustment does. The engine's t
matched statsmodels to within 1 × 10⁻¹⁴ in all three.

| target B (chr) | raw correlation r (naive t) | adjusted for B's own CN only | all covariates **except** arm CN | **full model t** | effect β | TRAIN hit? |
|---|---|---|---|---|---|---|
| MYBL2 (20) | 0.20 (7.8) | 8.4 | 7.3 | **7.49** | +0.65 | yes, down |
| NT5E (6) | −0.26 (−10.1) | −10.0 | −8.1 | **−7.35** | −1.18 | yes, up |
| GREB1 (2) | 0.16 (6.1) | 5.9 | 2.7 | **1.46** | +0.13 | no |

- **MYBL2:** more PPP1R15B copies go with *more* MYBL2 (a proliferation transcription factor). The association
  barely changes under adjustment.
- **NT5E** goes the other way, and its association shrinks somewhat once lineage and arm events are accounted
  for.
- **GREB1** is the warning case. Of all benchmark targets whose adjusted |t| is below 2, it is the one with the
  strongest raw correlation, so it was picked automatically. Its raw signal (t = 6.1) mostly disappears:
  - lineage and the other shared covariates take it to 2.7;
  - arm-level CN of other chromosomes takes it to 1.5.
  - A plausible reading: 1q gain is common in the lines (for example breast) that express this oestrogen
    target, and 1q is co-gained with other arms. That is cross-chromosome confounding, which the covariates are
    there to remove.

`step4_fwl_scatter` in the JSON has the raw and the fully adjusted (residualised) x_A vs y_MYBL2 points for a
before/after plot.

**Caveat (true of every 1q gene):** the model cannot separate PPP1R15B from the genes physically next to it,
because they gain and lose copies together. Across all 16,388 trans targets PPP1R15B's |t| has median 1.35,
about twice what pure noise would give (≈ 0.67). Much of that is a signature shared by the whole region. The
next step deals with this.

## Step 4 — Who in the co-amplified block gets the credit? (attribution A3)

**Segmentation.**
- Chromosome 1 has 1,938 genes with CN data.
- The tool walks along it and groups neighbours whose CN profiles (after covariate adjustment) correlate at
  r ≥ 0.95 with the running block average.
- This gives 483 blocks ("co-segments"). PPP1R15B's block has 4 genes within 0.3 Mb:

| gene | position (Mb) | own cis t | passes G1? | prior share π | median α across targets | perturbed in HCT116 screen? | TRAIN hits |
|---|---|---|---|---|---|---|---|
| **PPP1R15B** | 204.40 | **23.9** | yes | 0.854 | **0.856** | yes (TRAIN) | **835** |
| PIK3C2B | 204.42 | 11.6 | yes | 0.0018 | 0.0017 | yes (TRAIN) | 2 |
| MDM4 | 204.52 | 20.4 | yes | 0.144 | 0.142 | no | — |
| LRRN2 | 204.62 | 8.5 | **no** | 0.0004 | 0.0004 | no | — |

**How credit is shared** (a single-effect "SuSiE L = 1" posterior, computed separately for every target B). The
tool assumes exactly one gene in the block drives B and asks how likely each candidate is:

    α(g | B) ∝ prior(g) × evidence(g, B)
    prior(g)  ∝ exp(t_cis,g / 2)                       — genes whose CN clearly drives their own mRNA are favoured
    log evidence = ½·log(1−λ) + ½·λ·t_gB²,  λ = 0.05² / (0.05² + 1/df) ≈ 0.77

The arithmetic for B = MYBL2:

| gene | t with MYBL2 | log evidence | log prior | **α** |
|---|---|---|---|---|
| PPP1R15B | 7.49 | 20.87 | 11.96 | **0.934** |
| MDM4 | 7.34 | 19.99 | 10.18 | 0.066 |
| PIK3C2B | 7.05 | 18.41 | 5.80 | 0.0002 |
| LRRN2 | 7.14 | 18.86 | 4.26 | 0.0001 |

- All four genes correlate with MYBL2 about equally (t 7.0–7.5), because their copy numbers move together. The
  trans data alone cannot pick the driver.
- The decision therefore comes mostly from the **prior**: PPP1R15B's mRNA follows its CN most tightly
  (t_cis 23.9 vs 20.4 for MDM4).
- **Contrast neighbour 1, LRRN2:** it fails G1 (t_cis 8.5 < 10), so all its scores are exactly 0, and A3 also
  gives it almost no credit (α ≈ 0.0004).
- **Contrast neighbour 2, PIK3C2B:** it passes G1 (t_cis 11.6) and was itself knocked down in the screen, but A3
  gives it α ≈ 0.002. Its largest score is 0.08 (PPP1R15B's is 7.6), and none of its pairs reach the top 0.1 %.
  In the screen, PIK3C2B knockdown produced just **2** trans hits against PPP1R15B's 835. Here the attribution
  credited the right gene. This is one anecdote, not proof: attribution added +0.36 log-OR on TRAIN.
- **MDM4,** the well-known oncogene of this amplicon, keeps about 14 % of the credit. It was never knocked down
  in HCT116, so that share cannot be checked.

## Step 5 — Final score and predicted knockdown direction

    S(A, B) = t(A, B) × w_G1(A) × α(A | B)
    predicted change in B when A is knocked down = − sign(S) × s_A

For MYBL2: S = 7.49 × 1 × 0.934 = **7.00**. More copies go with more MYBL2, so knocking PPP1R15B down predicts
MYBL2 goes **down**. The screen saw MYBL2 go **down**, so the direction is concordant.

Across all of PPP1R15B's trans targets, α ranges from 0.70 (1st percentile) to 0.91 (90th percentile). So
PPP1R15B keeps most of the block's signal for every target.

### PPP1R15B's top 10 benchmark-testable targets, ranked by |S|

These are targets measured in the HCT116 screen.

| rank | target (chr) | t | α | S | predicted on knockdown | TRAIN hit? | observed | agrees? |
|---|---|---|---|---|---|---|---|---|
| 1 | POLD4 (11) | −7.67 | 0.99 | −7.61 | up | no | | |
| 2 | SLC39A13 (11) | −7.76 | 0.97 | −7.52 | up | no | | |
| 3 | MYBL2 (20) | 7.49 | 0.93 | 7.00 | down | **yes** | down | ✓ |
| 4 | NT5E (6) | −7.35 | 0.93 | −6.87 | up | **yes** | up | ✓ |
| 5 | GNAI2 (3) | −6.84 | 0.98 | −6.69 | up | no | | |
| 6 | ZNF273 (7) | 6.84 | 0.98 | 6.69 | down | no | | |
| 7 | MVP (16) | −7.24 | 0.92 | −6.67 | up | no | | |
| 8 | RRAS (19) | −6.94 | 0.96 | −6.66 | up | **yes** | down | ✗ |
| 9 | GSN (9) | −6.77 | 0.98 | −6.65 | up | no | | |
| 10 | ZNF711 (X) | 6.68 | 0.99 | 6.60 | down | no | | |

Over all trans targets, the top 10 also includes 4 genes that the benchmark cannot score (SPSB4, KLHL40, EVA1A,
MYBPC2). They are outside the HCT116 benchmark universe because the screen does not detect them consistently.
These are listed in `step6_top10_all_trans_targets`.

## Step 6 — How PPP1R15B fares in the benchmark

The benchmark ranks all 36.7 million TRAIN pairs by |S| and takes the top 0.1 % (K = 36,665 pairs, |S| ≥ 3.93).

| | value |
|---|---|
| PPP1R15B pairs in the top 0.1 % | **448**, 9th most of 3,806 TRAIN regulators (only 490 regulators have any) |
| hits among those 448 | 53 (**11.8 %**) |
| hit rate among PPP1R15B's other 8,731 testable targets | 9.0 % |
| hit rate across all TRAIN pairs (background) | **0.22 %** |
| within-PPP1R15B odds ratio (top-ranked vs other targets) | 1.36 |
| top-ranked hits whose predicted direction is right | 34 of 53 (64 %) |
| label-free rank by strongest score, all 12,141 regulators | 454th |

This one gene reproduces the project's main finding in miniature:
- **Choosing the right regulator does most of the work.** Pairs with PPP1R15B are hits about 40 times more
  often than random pairs (9–12 % vs 0.22 %). The G1 gate selects regulators like this one.
- **Ranking its targets helps only a little.** The top-ranked targets are hits only modestly more often than its
  other targets (odds ratio 1.36). Across all TEST regulators, within-regulator ranking showed no detectable
  information (AUROC 0.500, CMH 0.94).
- **Direction is right in 64 % of cases here, but whether that beats chance is not established.** On TEST, pooled
  concordance was 72 % against a 65 % within-regulator null (p = 0.012), but a permuted-CN negative control passes
  the same test (p = 0.001; correction D9), so the test is miscalibrated.

## Step 7 (optional) — The context layer

The reported, non-selected variant rescales scores towards genes that are highly expressed in HCT116, because
those are the changes Perturb-seq can detect:

    S_ctx = sign(t) × ρ_P × ρ_T(B)²

- ρ_P is the percentile of |S| among all 210 million trans scores. 55 % of them are exactly 0, from regulators
  that fail G1. MYBL2's |S| is exceeded by only 2,258 pairs.
- ρ_T is B's expression percentile in DepMap's HCT116 profile.

Under this layer, PPP1R15B's top 10 testable targets become ACTG1, TPT1, MT2A, RPS12, PTMA, EEF1A1, TMSB10,
ANXA2, RPS6 and RPL7. These are all among the 0.4 % most highly expressed genes in HCT116, although their
CN-only |S| is only 2.5–6.6.
- 7 of the 10 are TRAIN hits, against 3 of 10 for the CN-only ranking.
- 6 of those 7 have the right direction.
- What the layer adds is detectability: highly expressed genes are the ones Perturb-seq can measure. It adds
  little biology. On TEST the context layer reached OR 18.4, no better than a DepMap-only "hub" baseline (B5,
  22.8), and its pooled direction agreement was not significant (0.549 vs 0.507).

## What to take away

1. PPP1R15B's copy number varies naturally across cancer lines, mostly through 1q gains, and its mRNA tracks
   that dosage closely (t_cis = 23.9). That makes it a usable "natural experiment".
2. After adjusting for lineage, culture and arm-level co-gains, its CN still associates strongly with hundreds
   of genes elsewhere in the genome. Adjustment matters: GREB1 dropped from t 6.1 to 1.5.
3. Neighbouring genes on the same amplicon are statistically indistinguishable in trans. The attribution step
   splits credit using each gene's own dosage response, giving PPP1R15B ~85–93 % and MDM4 ~14 %.
4. In the knockdown screen, PPP1R15B is a strong trans regulator, and the tool ranks it near the top. Which of
   its targets respond, however, is predicted only slightly better than chance.
