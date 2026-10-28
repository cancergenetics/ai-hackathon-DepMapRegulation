# Turning points: how DepMapTR arrived at its final method

This page covers the decisions that produced the final method. It does not replay the whole research tree.
Times are UTC on 2026-10-27, taken from `git log --date=iso`. "Reviewer" means an AI reviewer agent; it never
means a human. Literature tags such as [F1] and [T1] are the reference keys used in `LITERATURE.md`.

**Final method.** The method has five steps and an evaluation design. The step IDs are used throughout these pages.

| id | step | final setting |
|---|---|---|
| S1 | association model | LOCO reduced form (leave-own-chromosome-out). Expression of B is regressed on the CN of A, adjusting for: B's cis-CN basis [x_B, homdel_B]; arm CN of every chromosome except A's; LOCO autosomal burden; lineage, ploidy/WGD, MSI, TP53, sex, media, growth pattern and RNA strandedness. Fitted by FWL on n = 1,436 lines (S2). |
| S2 | cis gate (G1) | Keep regulator A only if its own CN→mRNA first-stage F = t²_cis ≥ 100. Other regulators get weight 0. |
| S3 | attribution (A3) | SuSiE L=1 single-effect posterior over the greedy co-segment (r ≥ 0.95). Prior ∝ exp(t_cis/2), σ₀ = 0.05, no co-expression term. |
| S4 | scoring / sign | S = t·w_G1·α, ranked by \|S\|. Predicted knockdown sign = −sign(S)·s_A, where s_A = −1 only if t_cis ≤ −3. |
| S5 | context layer | Reported separately, not selected: S = sign(t)·ρ_P·ρ_T². ρ_T is B's expression percentile in DepMap HCT116. |
| E | evaluation design | Arm-level hold-out, sealed before tuning, with one LOCKed TEST run. Regulator-label permutation and arm jackknife. Hub and detectability baselines (B4e, B5). Within-regulator controls (CMH, C1, C2). Sign concordance against c₀ and a within-regulator null. |

---

## TP1 — Use the reduced form, adjust for the other chromosomes' arms, never for expression factors
- **When / where:** 18:21 to 18:33. Commits `499dc16` (LITERATURE.md, 59 refs) and `a17f5bd` (DESIGN.md v1,
  §4–§5.1). Refined in round 1 (`6d755e0`). Files: LITERATURE.md §2, §3, §6 and "Implications"; DESIGN.md §4, §5.1.
- **Situation / question:** The brief asks for a CN→trans-expression model that adjusts for cis CN and covariates.
  It also warns that correlated arm-level events can confound trans pairs. The question was how to remove that
  confounding without removing real trans signal.
- **Options considered:**
  - Regress out expression PCs or PEER factors, the eQTL default.
  - CN-derived PCs.
  - Arm CN of all arms except A's chromosome (LOCO).
  - Focal adjustment, which also removes A's own arm.
  - GLS with a CN kinship matrix.
  - A per-lineage meta-analysis.
- **Evidence that decided it (from the literature):**
  - Fehrmann 2015 [F1]: CN correlates positively with mRNA for 99% of abundantly expressed genes, so top
    expression PCs are partly CN effects.
  - Kang 2008 (ICE) [K1]: correcting for global expression correlation removes spurious trans hotspots.
    LITERATURE.md notes, citing "later work" without a specific reference, that it also removes genuine hotspots.
  - Dahl 2019 [D2]: adjusting for phenotype PCs causes >10-fold false-positive inflation, a collider-type bias.
  - Taylor 2018 [T1] and Shukla 2020 [S4]: recurrent arm events appear in >30% of tumours, against <15% for focal
    events, and there are 88 co-occurring arm pairs. Davoli [D3] adds that arms are co-selected, so cross-chromosome
    arm correlation is the main trans confounder.
  - "LOCO" is borrowed from GWAS practice. LITERATURE.md cites no specific paper for it.
  - Fast closed-form FWL / matrix-product partial correlations follow Matrix eQTL [S5] and tensorQTL [T2].
  - Reviewer-agent refinements in round 1:
    - burden computed leave-one-chromosome-out (domain #3);
    - TP53 and RNA strandedness added (domain #1–2, data #5);
    - cis basis with a homozygous-deletion indicator (statistics #2);
    - float64 QR plus NaN guards (algorithms #1–2).
- **Outcome:**
  - M1 is the LOCO reduced form, fitted with one pivoted-QR residualisation per regulator chromosome. It matches
    statsmodels OLS to 1e-8 on 480 real pairs.
  - Alternatives tested at tier 1 did not beat LOCO: focal +0.68 (unstable), leave-one-arm-out +0.04, CN-PCs(k=20)
    −0.02 (`results/m4/tier1_compare.txt`).
  - Expression factors were never used. A regulator-conditional hidden-factor sensitivity analysis, the GBAT [L2] /
    PANAMA [F2] "supervised factor" idea, scored TEST 2.90 against 3.85 (Δlog-OR −0.28, CI −0.64 to 0.07).
- **Caveat:**
  - The core premise was never tested directly: plain expression-PC adjustment was not run. "PCs remove real
    signal" is believed from the literature, not shown here.
  - M1 on its own was weak: TRAIN OR 1.59, permutation p 0.25.
  - On TEST, the paired gain of the selected method over naive Pearson (B1) is not significant (+1.22, CI −0.48 to
    2.92).
- **Produced:** S1.

## TP2 — Data profiling makes HCT116 the primary line and shows the benchmark is hub-dominated
- **When / where:** 18:32. Commit `aacc223` (docs/DATA_PROFILE.md §10–11). Decision D1 was recorded in the tree at
  `909a6b7` (18:34).
- **Situation / question:** Which of the two X-Atlas/Orion lines should be primary, and what does a positive
  actually look like?
- **Options:** HCT116, HEK293T, or both weighted equally.
- **Evidence:**
  - In HEK293T, 83.0% of trans hits come from perturbations whose own gene is not expressed; in HCT116 the figure is
    12.7%. In HCT116, non-expressed perturbations average 5.5 hits against 51.0, a plausible null level.
  - HEK293T is not a DepMap model; HCT116 is.
  - Hits are concentrated: the top 1% / 10% of perturbations give 55% / 92% of HCT116 trans hits.
  - LFC magnitude is a pseudocount artefact: in 83.5% of expressed targets the perturbed-cell mean is exactly 0.
- **Outcome:**
  - HCT116 became primary and HEK293T secondary (D1). This was decided before any model existed.
  - Labels use sign and p only; LFC magnitude is never used.
  - Non-expressed perturbations became an empirical null (class iii, NC2).
  - Hub concentration motivated the stratified metrics and hub baselines added in TP3.
  - Validation: on TEST the CN-only tool fails in HEK293T (OR 0.49, 6 overlaps), so making it secondary was right.
- **Produced:** E (primary line, label definition, null classes).

## TP3 — Round-1 review rebuilds the evaluation: arm hold-out, dependence-aware statistics, hub baselines, sealed TEST
- **When / where:**
  - 18:45: `6d755e0` (7 round-1 reviews, DESIGN v2, §13 dispositions).
  - 18:51: `9f7d21a` (round 2).
  - 18:55: `a6ccaeb` (v3.1 frozen).
  - 19:29: `57fbe32` (arm split built).
  - 20:05: `57f2de3` (TEST physically sealed).
  - Files: reviews/round1/*.md, DESIGN.md §6 and §13, results/split/.
- **Situation / question:** DESIGN v1 had four weaknesses:
  - 5 Mb-bin hold-out;
  - Fisher p-values that assume independent pairs;
  - variance-only hub baselines;
  - a PureCN-complete sample set (S1).
- **Options:**
  - Hold-out unit: 5 Mb bins, ≥25 Mb blocks, arms, or whole chromosomes.
  - Significance: Fisher, within-regulator score shuffle, regulator-label permutation, or a degree-preserving null.
  - Sample set: S1 or S2.
- **Evidence:**
  - Gene-gene CN correlation stays at r 0.56–0.82 at 5–10 Mb (data review #4). With bin splits, TEST regulators
    would have near-identical TRAIN neighbours.
  - Four reviewers proposed an arm split independently: AI-failure #4, external PI #6, prior work #4 and data #4.
  - At least 4 of 7 reviewers independently said Fisher p is invalid under pair dependence: external PI #2,
    statistics #9, prior work #2, AI-failure #10. The top 1% of regulators carry 55% of hits.
  - S1 silently dropped HCT116, the benchmark's own line, together with 45% of MSI-H lines and 31% of samples
    (data review E2).
- **Outcome:**
  - **Split.** The hold-out unit became the chromosome arm: 14 TEST and 28 TRAIN arms, assigned by triplets with
    seed 20261027. The leakage check: a TEST regulator's maximum CN correlation with any TRAIN regulator has median
    0.19, and none exceeds 0.8 (`results/split/cn_leakage_summary.json`).
  - **Significance.** The primary test became a regulator-label permutation (statistics #9a). It was preferred to the
    external PI's within-regulator shuffle because it preserves target hubness. CIs come from a delete-one-arm
    jackknife.
  - **Degree-stratified permutation.** A secondary permutation stratified by degree was added (round-2 statistics
    #3). It is a weaker version of the degree-preserving null recommended in LITERATURE §5 from Warmerdam 2026 [W1].
  - **NC1 null.** NC1 permutes the residualised CN matrix within lineage. Statistics reviewer #5 proposed this "as in
    tensorQTL". LITERATURE.md cites tensorQTL only for speed, so this idea came from the reviewer.
  - **Baselines.** B4e (target expression in HCT116; prior work #3) and B5 (combined hub; external PI #11) were added.
  - **Sample set.** S2 (n = 1,436) became primary.
  - **Sealing.** TEST labels and raw DE files were moved out of the tree and made unreadable, with their hashes
    committed (AI-failure #1–2, round-2 R1/R2).
- **Modified, not adopted:** Statistics #11 proposed tuning on the CMH OR. The PI kept the pooled OR because "the
  brief's primary is pooled enrichment". In hindsight this mattered (see TP6 and TP7).
- **Cost:** Only 14 TEST arms, so TEST precision is low. The final CI is 1.33–11.2.
- **Produced:** E, and the S1 sample set.

## TP4 — Reviews turn the PI's cis-instrument and attribution sketches into specified components
- **When / where:** 18:45 to 18:55. Commits `6d755e0` (v2), `9f7d21a` (v3) and `a6ccaeb` (v3.1). Files: DESIGN.md §1,
  §5.3, §5.4, §5.6; reviews/round1/statistics.md #6–8; reviews/round1/algorithms_compute.md #8–9;
  reviews/round3/external_pi_fresh.md #1–3.
- **Situation / question:** v1 sketched three ideas:
  - an instrument gate, "M2", with a tuned τ;
  - a co-expression "mediation" score;
  - a SuSiE-like softmax over co-segmented genes.

  The predicted knockdown sign was simply sign(β_AB).
- **Options:**
  - Gating: a gate or a 2SLS score.
  - Attribution: mediation, a softmax, iterative SuSiE (IBSS), or closed-form single-effect fine-mapping.
  - Sign: the raw reduced-form sign, the IV sign, or a shrunk cis sign.
- **Evidence and credit:**
  - **Instrument logic.**
    - From the literature: Goh's "genomic randomisation" / Mendelian-randomisation framing [G2] and weak-instrument
      bias (Burgess & Thompson 2011 [B2]).
    - From the literature: the only strong published enrichment is cis-anchored, eQTLGen2 vs X-Atlas HCT116 OR 12.4
      (Warmerdam 2026 [W1]). That figure rests on 41 pairs with a degree-matched null.
    - Reviewer: prior work #1 asked for the gated score to be first in line.
    - Reviewer: statistics #6 pointed out that the reduced-form t is the Anderson–Rubin test, and proposed an F ≥ 10
      floor with τ tuned at or above it. The {10, 30, 100} grid itself first appears in DESIGN v2, written by the PI.
  - **Attribution.**
    - SuSiE [W4] comes from LITERATURE.md.
    - The **closed-form Wakefield ABF with L = 1** came from the algorithms/compute reviewer (#9); it is not in
      LITERATURE.md. That reviewer also warned that CN-only z cannot discriminate co-segmented genes at r ≈ 0.97–1.0,
      so discrimination must come from the prior.
    - Statistics #8 also proposed L = 1.
    - Round 2 fixed σ₀ = 0.05 a priori (statistics #5) and normalised over all segment genes so the benchmark could
      not leak in (external PI #3).
  - **Mediation demoted.** Co-expression "mediation" was renamed residual co-expression and confined to within-segment
    use (statistics #7, prior work #6, citing Aguirre 2025 [A1]: across pairs, co-expression is uncorrelated with
    perturbation effects).
  - **Sign.**
    - The IV sign came from domain #7 and statistics #6.
    - Round 2 shrank it: s_A = −1 only if t_cis ≤ −3, because about 95% of genes scale positively (domain #1).
    - Round 3 caught that §1 and §6 contradicted each other on the sign rule and that the final score was undefined.
      It defined S = t·w_G1·w_G2·α (round-3 #1, #3).
- **Outcome:** G1, G2 and A3 became pre-registered, hash-locked variants with fixed grids (`results/variants.yaml`).
  The sign rule was unit-tested, including the case of a negative-cis regulator.
- **Produced:** S2, S3, S4.

## TP5 — First TRAIN numbers: a DepMap-only hub prior beats the CN model 12-fold
- **When / where:** 20:08. Commit `f2deb85` (`results/m4/mvp_full_hct116.txt`, `mvp_quick.txt`).
- **Situation / question:** Does the frozen core model (M1) enrich Perturb-seq pairs at all?
- **Evidence (HCT116 TRAIN, OR at top 0.1%):**

  | score | OR | note |
  |---|---:|---|
  | M1 | 1.59 | perm p 0.25 |
  | B5 (regulator cis×CN-SD × target HCT116 expression) | 19.7 | |
  | B4e (target expression) | 7.1 | |
  | B2 (co-expression) | 5.8 | |
  | B0 random | 0.98 | control |
  | NC1 | 0.75 | control |
  | B6 arm-only | 0.97 | control |

  - The frozen CMH stratification does not remove context detectability: B4e, with no pair information, scores
    CMH 3.8.
  - M1 sign concordance is 0.581 against c₀ 0.476, but the within-regulator sign permutation gives p = 0.66.
- **Outcome:**
  - The PI did not stop or redesign. It ran the pre-registered tiers as written and recorded that most pooled
    enrichment comes from hubness and detectability.
  - LESSONS #1 records the cost: B5 was defined at round 1, but its number only appeared after the design was frozen
    and sealed.
- **Produced:** E (interpretation framework). This result directly motivated the TP7 controls.

## TP6 — The cis-instrument gate is the one large, stable gain; A3 helps through its cis prior
- **When / where:** 20:21 (`6288007`, tiers 1–2), 20:37 (`c968b1d`, tiers 3–4) and 20:47 (`48143dc`, tiers 5–7).
  Files: results/tuning/tier2_G1_compare.json, tier3_G2_compare.json, tier4_A3_compare.json; RESEARCH_TREE M5.
- **Situation / question:** Which pre-registered component survives the strict adoption rule (DESIGN §8.2: both
  folds > 0 and paired 90% lower bound > 0)?
- **Evidence (HCT116 TRAIN, Δlog-OR vs incumbent):**
  - **G1:** τ = 10 gives +0.25, τ = 30 +0.51, τ = 100 +0.79 (folds 0.65 / 1.12, jackknife lb90 +0.61). The gain rises
    monotonically up to the grid maximum. 5,403 of 12,141 regulators pass. ΔCMH −0.05, so the gain is regulator
    selection, not target ranking. With τ = 100, w = 0 and w = 0.5 give identical results: gated regulators fill
    the whole top-K.
  - **G2** (A expressed in HCT116) on top of G1: +0.05, lb90 ≈ 0. Not adopted; regulators that pass a strong cis gate
    are already expressed.
  - **A3:** ρ 0.95 / s 2 / w_co 0 gives +0.36 (folds 0.43 / 0.28, lb90 +0.04). w_co 0.25 also passes (+0.34);
    ρ 0.98 (−0.05) and s 5 (0.00) do not.
  - **Tiers 5–7** all fail: target standardisation −0.03, hinge −0.11, floor 2⁻⁵ −0.01, RINT −0.16.
  - **Incumbent:** TRAIN OR 5.0 (perm p at its 1/201 floor), CMH 1.09.
- **Outcome:**
  - G1 and A3 were adopted. LOCO was retained.
  - Selection used 23 of the 30-configuration cap.
  - The fresh-scientist diagnostic (TP7) later showed that A3 is nearly a hard argmax. Its gain is regulator
    selection through the exp(t_cis/2) prior, not attribution of targets (α within-regulator AUROC 0.498).
- **Produced:** S2 (τ = 100, w = 0), S3 (ρ = 0.95, s = 2, w_co = 0) and confirmation of S1 (LOCO).

## TP7 — Diagnosis: CN carries no within-regulator target information, so the detectability fusion becomes a separate context layer
- **When / where:**
  - 20:47: `d4bdd91` (retro).
  - 21:02: `704aac0` (D6; reviews/ideas/M5a_ideas.md).
  - 21:48: `f209390` (I1–I6, C1/C2 controls).
  - 22:02: `d439c58` (D6a loader fix, D7 selection).
  - Files: results/tuning/d6_*.json; RESEARCH_TREE D6, D6a, D7.
- **Situation / question:** The incumbent sits far below B5 (Δlog-OR −1.37, CI −1.91 to −0.83). Should the remaining
  7 configurations chase pooled OR, and how should any gain be credited?
- **Evidence (from a fresh "ideas scientist" agent, using TRAIN only):**
  - Within a regulator, \|t\| ranks targets at noise level: AUROC 0.507. Target expression in HCT116 reaches
    0.70–0.71 and beats \|t\| for 96% of regulators.
  - About 94% of the log-enrichment is regulator allocation; the pair term contributes about 1.09×.
  - A G1 regulator's t-row correlates with its arm-only row at median r = 0.87.
  - There is a negative-t component that scales with target expression. It is partly a TPM compositional artefact.
  - Proposed controls:
    - C1: shuffle each regulator's row within 50 HCT116-expression bins;
    - C2: apply the same fusion to permuted-CN t.
- **Options tried (D6, post hoc, label-free):**
  - I5 compositional renormalisation: +0.05;
  - I4 coarser segments: +0.21, unstable;
  - I3 shrunk-IV score: +0.18, ΔCMH +0.42;
  - I6 lineage replication: +0.04;
  - I1 detectability fusion: γ = 1 gives +1.50, γ = 2 gives +1.57 (folds 1.37 / 1.82), TRAIN OR 24.1.
- **Decisive evidence:**
  - I1 passed §8.2, but the exact C1 gave OR 24.5, so Δlog-OR(pair) = −0.02 (CI −0.09 to 0.06).
  - C2 (fusion on permuted CN) gave 13.9 against 24.1, so the real CN evidence does improve regulator allocation.
  - Against B5: +0.20 (folds 0.04 / 0.41, lb90 < 0), not significant.
  - D6a: before a loader fix, C1 ran in a biased approximate mode and reported Δ +0.063 (lb90 +0.002), a spurious
    "pair-specific gain".
- **Outcome (D7):**
  - The selected method stays CN-only (S1–S4).
  - The I1 γ = 2 fusion ships as a separately named **HCT116 context layer**, so the estimand output is not mixed
    with a detectability prior.
  - I2 was not run, because its adoption rule could not be satisfied by construction.
- **Credit and missed chance:**
  - The fusion and the C1/C2 controls came from the ideas-scientist agent.
  - Round-1 prior-work reviewer #5 had already proposed this target-expression prior and called it "likely the
    largest OR gain" (it is "hub exploitation").
  - DESIGN §13 dispositioned that proposal as "covered by B3/G1". That disposition is inaccurate, because B3 and G1
    cover only the regulator side. Had it been pre-registered, the context layer would not be post hoc.
- **Produced:** S5, and E (C1/C2 within-regulator controls; the HCT116-stratified CMH diagnostic).

## TP8 — One sealed TEST run; the cis gate replicates, the headline is not significant, and the claims are narrowed
- **When / where:**
  - 22:54: `7b4dba9` (pre-final review; reviews/code/prefinal_review.md).
  - 23:07: `5d92561` (D8a fixes).
  - 23:08: `62e67b9` (selected.yaml), then `603a7f5` (LOCK STARTED).
  - 23:12: `e9085c7` (TEST COMPLETED).
  - 23:26: `b0827fd` (REPORT corrected after verification).
  - 23:28: `c27a904`.
  - Files: results/final/summary.md; reviews/report/report_verification.md.
- **Situation / question:** Does the TRAIN selection hold on held-out arms, and what may be claimed?
- **Pre-run safeguards:** A fresh pre-final reviewer gave a conditional GO with three fixes:
  1. sealed-file hash preflight;
  2. a NaN CMH bootstrap bound;
  3. the missing context-layer comparisons.

  A mock run on TRAIN arms relabelled as TEST exercised the same code path.
- **Evidence (HCT116 TEST: 14 arms, 1,817 regulators, 17.7 M pairs):**
  - **Selected tool:** OR 3.85 (jackknife CI 1.33–11.2). Regulator-label permutation p = 0.068, which is **not
    significant**. Gain over NC1 +1.10 (CI −0.66 to 2.85).
  - **G1 ablation:** 1.46 → 3.85, Δlog-OR +0.97 (CI 0.30–1.64). This replicates TRAIN (+0.91 for the same removal
    contrast; +0.79 when G1 was added to M1).
  - **Within regulator:** AUROC 0.500, CMH 0.94 (CI 0.77–1.15), C1 3.76.
  - **Sign:** 0.718 against a within-regulator null mean of 0.647, p = 0.012, n = 131. It is not significant in S1
    (p 0.13) or in the no-split-sibling subset (p 0.48).
  - **Context layer:** 18.4 (CI 10.9–30.9), against B5 22.8 (Δ −0.22, CI −2.3 to 1.9) and C2 10.6 (Δ +0.55,
    CI 0.07–1.03).
  - **Report verification:** a fresh verifier checked about 120 numbers and found 4 errors and 3 over-claims, all
    fixed.
- **Outcome:**
  - The final claims are narrowed to: CN may help identify **which regulators** have detectable trans consequences
    (resting mainly on the replicated G1 gain); it showed no detectable information about **which targets** respond.
    A "modest" direction claim made at the time was withdrawn under D9 (28 Oct): the C2 control passes the same
    sign test (p = 0.001), so direction is not established.
  - Precision is limited by 14 arms.
- **Produced:** E (one-shot protocol and reporting rules) and the evidential status of S2 and S4.

---

### Appendix: who contributed which idea to the final method

| idea | first source | type |
|---|---|---|
| No expression PCs/PEER | Fehrmann [F1], Kang [K1], Dahl [D2] | literature |
| LOCO arm CN | Taylor [T1], Shukla [S4], Davoli [D3] (motivation); "LOCO" from GWAS practice (no citation); the brief's own warning | literature + human brief |
| LOCO genome burden; TP53; RNA strandedness | domain reviewer R1 #1–3 | reviewer agent |
| Cis basis with homdel indicator | statistics reviewer R1 #2 | reviewer agent |
| Cis-instrument gate (G1) | Goh [G2] MR framing, Burgess [B2], Warmerdam [W1]; PI v1 "M2"; prior-work R1 #1; statistics R1 #6 (F, AR) | literature + PI + reviewers |
| SuSiE L=1 attribution | SuSiE [W4] (literature); L=1 + closed-form **Wakefield ABF**: algorithms R1 #9, statistics R1 #8 | literature + reviewer (ABF not in LITERATURE.md) |
| Shrunk cis sign (t_cis ≤ −3) | domain R2 #1, statistics R1/R2 #6 | reviewer agent |
| Arm-level hold-out | AI-failure #4, external PI #6, prior work #4, data #4 (R1) | reviewer agents (supported by [T1]) |
| Regulator-label permutation | statistics R1 #9a | reviewer agent |
| Degree-stratified permutation (secondary) | statistics R2 #3; LITERATURE recommends Warmerdam-style degree-preserving null [W1] | reviewer, partly literature |
| NC1 residual permutation | statistics R1 #5 ("as in tensorQTL") | reviewer agent (tensorQTL [T2] cited in LITERATURE only for speed) |
| B4e / B5 baselines | prior work R1 #3 / external PI R1 #11 | reviewer agents |
| Physical TEST seal + LOCK | AI-failure R1 #1–2, R2 R1/R2 | reviewer agent |
| C1 / C2 controls; context fusion I1 | ideas scientist (M5a); target prior first proposed by prior work R1 #5 | agents |
| Regulator-conditional hidden factors (sensitivity) | GBAT [L2] / PANAMA [F2]; prior work R1 #7 | literature + reviewer |
