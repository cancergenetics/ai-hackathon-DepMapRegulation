# M5a method ideas: fresh-scientist review (2026-10-27)

Scope: HCT116 TRAIN only. The incumbent is `g1_t100_w0__a3_r095_s2_w0` (LOCO + G1 τ=100 + A3 ρ=0.95, s=2). Inputs read:
PROMPT, DESIGN v3.1 (§8), RESEARCH_TREE (M5a), docs/MODEL.md, docs/BENCHMARK.md, results/tuning/*, results/eval/*.

**Label hygiene.**
- No `evaluate()` or `compare()` calls were made, so nothing was appended to either tuning log.
- Label-touching diagnostics were direct computations with `benchmark.Evaluation` on TRAIN positives (diag1–5, diag8). diag6 and diag7 are label-free.
- No candidate configuration was scored or tuned. The one exception is the decomposition bound in §1.4, flagged there.
- Scripts are in `reviews/ideas/M5a_diag/diag*.py`. Each runs in under 10 s with `.venv/bin/python`.

---

## 1. Diagnosis: why the CN pair evidence barely ranks targets within a regulator

### 1.1 Within-regulator ranking by |t| is at noise level, and target expression dominates it
Mean within-regulator AUROC over TRAIN regulators with ≥ 10 positives (n = 333; 228 pass G1), from diag2, diag3 and diag8:

| ranking of a regulator's targets | AUROC |
|---|---|
| incumbent \|t\| (LOCO) | **0.507** (SE ≈ 0.004); G1 regulators 0.510 |
| focal \|t\| (own arm adjusted) | 0.511 |
| \|S\| including A3 α | 0.508; α alone 0.498 |
| B2 co-expression \|r\| | 0.502 |
| B6 arm-only \|t\| | 0.508 |
| CN target responsiveness (mean \|t\| over regulators) | 0.532–0.536 |
| \|t\|·SD(y_B) (≈ effect size) | 0.540 |
| **target expression in DepMap HCT116 (B4e)** | **0.700–0.709** (beats \|t\| for 96% of regulators) |
| eh + 0.25·rank\|t\| | 0.701 (adding \|t\| *lowers* AUROC) |
| \|t\| ranked within HCT116-expression deciles | 0.506 |

This does not depend on the regulator:
- by degree: AUROC 0.496 / 0.525 / 0.512 / 0.505 across degree bins [10,30), [30,100), [100,300), ≥ 300;
- by instrument strength: ρ(AUROC, cis t) = 0.02.

### 1.2 The top-K is set by regulator selection; the pair term adds about 8%
From diag1 and diag4:
- Top-K precision is 1.094%, against a base rate of 0.222%.
- If each selected regulator's top-K targets were random, precision would be 0.985%. If they were random within HCT116-expression deciles, it would be 1.008%.
- So about 94% of the log-enrichment is regulator allocation. The pair-specific lift is about 1.09×, which matches CMH ≈ 1.09–1.10.

### 1.3 Top-K slots go by arm-signature magnitude, not by regulator relevance
From diag1, diag2 and diag4:
- **Arm concentration.** 490 regulators fill K = 36,665, and 100 of them take 71%. Arms 1q and 1p take 31% of K, yet their pair hit rates (0.0017 and 0.0025) are at base. 13q has a hit rate of 0.0057 but only 5% of K.
- **Rows are the arm signature.** A G1 regulator's t row correlates with its B6 arm-only row at median r = **0.87**, and 91% of rows exceed 0.7. B6 alone has OR 0.97.
- **Same-arm redundancy.** The 54 1q regulators in the top-K share only 1,308 distinct targets (redundancy 4.9).
- **Non-expressed regulators.** G1 measures the instrument across all of DepMap, not in HCT116:
  - MYCN (HCT116 log2(TPM+1) = 0.05, own class iii) takes 1,701 top-K slots for 1 hit;
  - MYCL (1.8) takes 618 slots for 10 hits;
  - KTI12 (class iii) takes 646 slots for 0 hits.
  - G1 regulators with context expression < 1 hold 8.2% of K, and those < 2 hold 13%.
- **Steep gradient in regulator context expression.** The per-pair hit rate of G1 regulators rises with HCT116 expression:

  | HCT116 log2(TPM+1) of A | < 1 | 1–3 | 3–5 | 5–7 | ≥ 7 |
  |---|---|---|---|---|---|
  | hit rate per pair | 0.0006 | 0.0012 | 0.0036 | 0.0074 | 0.0124 |

  G2 tested only hard thresholds at 1 and 2, which catch just the first column.
- **Correlates of regulator degree.** Benchmark degree correlates with cis t (ρ = 0.23) and context expression (0.20 within G1). It does not correlate with CN variability (n lines with |log2CN| ≥ 1: ρ = −0.01). Yet |t| grows with CN variance, so the allocation rewards frequently altered arms (1q gain, MYCN amplicons) for no benchmark benefit.

### 1.4 Top targets ignore detectability
From diag1 and diag4:
- **Top-K targets look like the universe.** Their median HCT116 expression is 4.67, against 4.74 for all universe targets.
- **Positives are skewed high.** Their median is 6.01. The hit rate per target is about 4 for targets at 2–4 log2(TPM+1) and 23 for targets ≥ 7.
- **Detectability is target-level, not pair-level.** ρ(target degree, HCT116 expression) = 0.45, while ρ(target top-K count, target degree) = 0.02.
- **Decomposition bound (label-free procedure, evaluated once on TRAIN labels).** Keep the incumbent's per-regulator top-K counts but take each regulator's targets by HCT116 expression. Precision rises from 1.09% to **4.97%**, about **OR 24**, compared with B5 at 19.7. Treat this as an upper-range sanity number, not as a tuned result.

### 1.5 A systematic negative-t component that scales with target expression (label-free)
From diag6 and diag7:
- **Skew.** Among trans pairs of G1 regulators, t ≤ −3 is twice as common as t ≥ 3 (3.5% vs 1.7%). 69% of rows have a negative mean.
- **Scaling with expression.** The mean t per target falls monotonically with target expression quintile: −0.075, −0.155, −0.247, −0.291, −0.311.
- **Link to transcript share.** The negative row mean is stronger when x_A tracks the TPM share of A's own chromosome: ρ = −0.48 over 1,570 regulators on chromosomes 1, 8, 12, 13 and 20.
- **This is partly a TPM compositional artefact.** Re-normalising y to exclude the transcripts of chr(A) moves the chr1 row-mean t from −0.12 to +0.09, chr8 from −0.59 to −0.34, and chr13 from −0.67 to −0.57.
- **The rest is a shared own-chromosome aneuploidy response.** The LOCO burden excludes chr(A) by design.
- **Little effect on ranking.** The renormalised rows correlate with the originals at r ≥ 0.99.
- **Top-K pairs.** Negative-t pairs have lower precision than positive-t pairs: 0.95% vs 1.43%, within-regulator MH OR 1.42 (diag5). Across all |t| ≥ 3 pairs the hit rates are equal (0.0066 vs 0.0069), so this is a tail effect.

### 1.6 Synthesis
For a given A, the LOCO reduced-form t row is mostly two things:
- a signature of A's arm or chromosome, shared by its neighbours;
- a generic component that grows with target expression, part compositional and part aneuploidy response.

Neither is specific to A, and focal adjustment does not uncover a hidden gene-specific part (AUROC 0.511).

The HCT116 hits are hub targets set mainly by detectability and by generic knockdown responses. In this benchmark, then:
- pooled OR is earned through regulator selection (G1, A3's cis prior) and could be earned through target detectability;
- pair-specific knowledge is close to zero.

One thing does carry pair information. The *sign* of t: concordance is 0.561 vs c₀ 0.483, with within-regulator permutation p = 0.005. The CN evidence knows something about direction even though it cannot pick targets. The primary metric is direction-agnostic, so this cannot be used for it.

---

## 2. Proposals (DepMap-internal, label-free; ≤ 2 configurations each; total 7)

Ranking criterion: expected pooled-OR gain first. Pair-specific gain is given separately. "Pair-specific" means the gain survives the matched controls C1 and C2 below, not the CMH, for the reason in §3.1.

**Matched controls (evaluation-only, not candidate configurations; log them as controls with a decision note).**
- **C1 (pair shuffle).** For each regulator A, permute its score row across targets within 50 bins of HCT116-expression percentile (seeded, `SeedSequence(20261027)`), then apply the same fusion.
  - It preserves regulator allocation and target detectability and destroys only pair evidence.
  - "Pair-specific gain" = OR(fused) − OR(fused C1), with a paired jackknife.
- **C2 (NC1 fusion).** Run the same fusion on the NC1 (index 0) t matrix, keeping G1 and setting α = 1. This tests the CN evidence at both the pair and the regulator level, beyond the G1 gate.

### I1. Context-detectability fusion (target side), 2 configurations
- **Definition:**
  - ρ_P(A,B) = average-rank percentile of |S_inc,AB| among all finite trans pairs of the full tool matrix (12,141 × 18,276). It is 0 where S_inc = 0, i.e. where G1 fails.
  - ρ_T(B) = percentile of DepMap HCT116 (ACH-000971, the `context` parameter, as for G2 and B4e) log2(TPM+1) among **all 18,276 DepMap measured targets**. Do not rank within the PIE universe, because its membership is benchmark-derived.
  - **S = sign(t_AB) · ρ_P · ρ_T^γ, with γ ∈ {1, 2}**, and s_cis unchanged.
- **Why it should work:**
  - Within a regulator, target expression has AUROC 0.70 and |t| has 0.51.
  - The decomposition bound (§1.4) suggests OR around 20–24, about +1.4 to +1.6 log-OR over 5.0.
  - Unlike B5, the regulator side keeps the incumbent's G1/A3 allocation.
  - Physically, Perturb-seq power for a fixed LFC grows with the target's count depth, so this is a prior on *detectability of the predicted consequence* in HCT116.
- **Pair-specific gain:** expected ≈ 0, since adding |t| to eh *lowers* within-regulator AUROC. The CMH will rise mechanically (see §3.1). The real test is C1.
- **Drift risk: high.** The fused score answers "which CN-predicted consequences are detectable in HCT116", not β_AB. Ship it as a separately named context layer. The CN-only score stays the estimand output. Make "selected minus I1" one of the three TEST ablations.

### I2. Graded regulator-context prior (continuous G2), 1 configuration
- **Definition:** S ← S · ρ_R(A), where ρ_R(A) = percentile of HCT116 log2(TPM+1) of A among the 12,141 tool regulators. CN-only regulators get ρ_R = 0. Apply it on top of the I1 winner, or on the incumbent if I1 is not adopted.
- **Why:**
  - The per-pair hit rate rises about 20× across context-expression bins inside G1 (§1.3).
  - 17.5% of K sits on regulators with expression < 3, whose per-pair hit rates are 0.06–0.12%. For regulators ≥ 5 the rates are 0.74–1.24%.
  - The hard G2 (thresholds 1 and 2) only removes the bottom bin, which explains its +0.05.
  - The rationale is the one already accepted for G2: a knockdown or dosage change of a weakly expressed gene changes little in that context.
- **Expected:** pooled +0.2 to +0.5 log-OR. Pair-specific: 0.
- **Drift: moderate.** It is a context-relevance weight. Report it separately, as for G2.

### I3. Power-aware effect-size score (CN side), 1 configuration
- **Definition:**
  - Replace t with a shrunk IV effect: β̂_IV,AB = r_AB·σ̃_B / (r_cis,A·σ̃_yA). Here σ̃_B is the residual SD of y_B after D_c and W_B, and σ̃_yA is the residual SD of y_A after D_c, both available from the group residualiser. The ratio is well conditioned because G1 enforces F ≥ 100.
  - Shrink with a positive-part James–Stein factor: β̃ = β̂_IV · max(0, 1 − 4/t_AB²).
  - Score S = β̃ · w_G1 · α, ranked by |S|. Use it in ρ_P of I1, or alone if I1 is not adopted.
- **Why:**
  - |t| ∝ |r|·√df rewards regulators whose CN varies a lot (1q, MYCN), and CN variability does not predict degree (ρ = −0.01).
  - The knockdown estimand is the effect *per unit of A's mRNA* (β_IV × Δlog A), which removes that inflation.
  - Within a regulator, |t|·SD_B beats |t| (AUROC 0.540 vs 0.510).
  - It is the DESIGN §5.3 IV quantity, so it is estimand-faithful.
- **Expected:** pooled ±0.3 (it reallocates regulators; the sign is uncertain). Pair-specific: small +.
- **Drift: low.** It is the CN→expression estimand, rescaled.

### I4. Coarser co-segment attribution, 1 configuration (optionally 2)
- **Definition:** A3 with ρ = 0.90 (optional second configuration ρ = 0.85); s = 2 and w_co = 0 unchanged.
- **Why:**
  - The A3 sweep moved in this direction (0.98 < 0.95).
  - Arm-mate redundancy is large (1q: 54 regulators share 1,308 targets).
  - With σ₀ = 0.05 and s = 2, attribution is nearly hard (§3.4). Larger segments therefore give one cis-prioritised owner per target per segment, which frees top-K slots for other arms.
- **Expected:** pooled +0.1 to +0.3. Pair-specific: ≈ 0 (α has within-regulator AUROC 0.498).
- **Drift: low** for the estimand. The attribution accuracy of individual genes may fall: report the SIM recovery at ρ = 0.90 first. It is label-free, so it costs no budget.

### I5. LOCO compositional renormalisation of expression, 1 configuration
- **Definition:** in group c, replace y_B by log2(TPM_B / (1 − share_c) + 1), where TPM = 2^y − 1 and share_c = Σ_{g∈chr c} TPM_g / Σ_g TPM_g, per line, over protein-coding columns. It has no df cost, and the same transform applies to the cis model and to B6. ~20 lines in `pipeline.prepare`/`GroupCache.get`, since Y depends on the group.
- **Why:**
  - It removes a technical path, CN(A) → transcript share of chr(A) → TPM_B, that the LOCO model attributes to A (§1.5).
  - That path biases t negatively and most strongly for highly expressed targets, which are where the benchmark hubs are.
- **Expected:** pooled ±0.1 (rows correlate at r ≥ 0.99). Pair-specific: small + (better sign behaviour).
- **Drift: negative.** It *improves* estimand fidelity. Recommended even if it is neutral on OR. It needs a decision note, because it changes y.

### I6. Lineage-replicated (conjunction) pair evidence, 1 configuration
- **Definition:**
  - Split the 32 Z lineage levels into two disjoint halves H1 and H2 with a seeded greedy balance on n. Fix the split before any evaluation and commit it.
  - Fit M1 in each half; drop Z dummies that are constant within the half.
  - t_rep = sign(t₁)·min(|t₁|, |t₂|)·1[sign t₁ = sign t₂].
  - S = t_rep · w_G1 · α, with α and G1 from the full sample.
- **Why:**
  - This is the only proposal aimed at *pair* evidence. It treats replication across lineage contexts as a filter for confounded or lineage-specific trans associations (CRC 18q−/20q+, the sex and lineage programmes).
  - A dosage effect mediated by A's product should generalise.
- **Expected:** pooled −0.2 to +0.1 (it loses power). Pair-specific: uncertain, small +. Compute: 2 engine runs, about 2 × 70 s.
- **Drift: low.** Same estimand, stricter evidence.

### Suggested use of the 7 remaining configurations
- **Steps 1–3, CN core first (estimand-faithful), greedy on the incumbent:**
  1. I5 (1 configuration);
  2. I4 ρ = 0.90 (1 configuration);
  3. I3 or I6 (1 configuration each; run both if budget allows).
- **Step 4, context layer last:**
  - I1 γ ∈ {1, 2} (2 configurations), then I2 (1 configuration);
  - each judged by §8.2 **plus** C1 and C2;
  - report both the CN-only and the fused selection.

**Not recommended:**
- focal under G1: focal t has AUROC 0.511, and only the regulator allocation changes;
- G1 τ = 300: the gain is regulator-hub selection, which I2 does more directly;
- CN-derived target-responsiveness priors: 0.53, below eh;
- any co-expression weight;
- signed asymmetry weights: the effect is only in the tail and not stable.

---

## 3. Things that look wrong or misleading in the current pipeline

1. **The CMH is not net of HCT116 detectability.**
   - **Problem.**
     - The strata are regulator × deciles of *S2-mean* expression.
     - B4e, which has no pair or regulator information, scores CMH 3.77, and B5 scores 2.77.
     - With HCT116-expression deciles (diag2), B4e still scores 3.38 and B5 2.57, because the top decile is wide.
   - **Consequences.**
     - The requirement "beat B5 on CMH" and the §8.2 CMH guard can both be met by a target-expression prior alone. Any I1-type fusion will pass them trivially.
     - B5's CMH of 2.77 is not evidence that B5 knows anything pair-specific: its degree-stratified permutation gives p = 0.27 and OR/null 1.17.
   - **What to report instead.** By the degree-stratified permutation, the incumbent (OR/null 1.47, p = 0.055) already has *more* pair-specific signal than B5. Adopt C1, the degree-permutation ratio and the expression-adjusted within-regulator AUROC as the pair-specific evidence. This needs a decision note; the frozen metric is unchanged and these are added diagnostics.
2. **Compositional TPM artefact in the LOCO t** (§1.5). It is a known limitation (round-1 domain #13), but it is not small: 2:1 negative skew, scaling with expression. I5 fixes the technical part. The remaining own-chromosome aneuploidy response should be stated as part of the "A's arm dosage" caveat in DESIGN §1.
3. **G1 is a population instrument, not a context instrument.**
   - MYCN and MYCL pass with high cis t from neuroblastoma and SCLC amplicons, and are not expressed in HCT116.
   - They consume about 6% of K for 11 hits.
   - G2 was evaluated before A3 and at thresholds where only 53 G1 regulators fall. Revisit it as I2.
4. **A3 is nearly a hard argmax.**
   - With σ₀ = 0.05 and df ≈ 1,335, λ ≈ 0.77, so log ABF ≈ 0.385·z². Inside a segment with CN r ≥ 0.95, a Δz of about 1 at z ≈ 8 gives a 20× evidence ratio, and the cis prior exp(t_cis/2) is also steep.
   - 74% of G1 trans pairs have α < 0.5; every top-K pair has α ≥ 0.5.
   - α carries no within-regulator information (AUROC 0.498). The +0.36 gain is regulator selection through the cis prior, not attribution of targets.
   - This is not a bug in the formula, but the description "attribution" overstates what it does. Say so in REPORT.
5. **Headline permutation p values use n_perm = 200** (incumbent_full.json, B5). "perm p = 0.005" in RESEARCH_TREE is the floor 1/201, not an estimate. The final evaluation uses 1,000, which is fine; reword the tree.
6. **B3 is degenerate under top-K.** It is target-agnostic, so K covers whole rows of about 4 regulators: fold-2 log-OR −3.78 and OR@0.01% = 0.06. B3 is a poor yardstick for "regulator selection"; use B5's degree-null instead.
7. **Fusion percentiles must not be computed within the PIE universe** (see I1). The target universe is defined by benchmark features with control mean ≥ 0.05. Use DepMap-wide percentiles so that the tool stays label-free.
8. **Interpretation risk for REPORT.**
   - As currently specified, the best achievable pooled OR on this benchmark comes almost entirely from DepMap-only regulator and target priors.
   - Any headline number above about 10 should be shown side by side with its C1 control and with the CN-only score.
   - Otherwise the result reads as "CN predicts Perturb-seq pairs" when it is "DepMap context predicts Perturb-seq hubs".
