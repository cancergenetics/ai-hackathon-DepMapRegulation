# Round 1 review: statistics (model and inference)

Reviewer remit: DESIGN.md v1, §3–§6. Notation: tilde = residual after projecting on D_c = [1, Z, arms∉c], for regulator chromosome c.

## Verdict on the core FWL scheme

The scheme is **algebraically correct**. X and Y are both residualised on the same D_c, and B's arm is in arms∉c because chr(B) ≠ c. So x̃_A, ỹ_B and x̃_B are residuals from one projection. The recursion
r_{AB·B} = (r_{AB} − r_{A,xB} r_{B,xB}) / √((1−r²_{A,xB})(1−r²_{B,xB}))
then conditions on exactly Z ∪ arms∉c ∪ {x_B}, as claimed. Collinearity between x_B and arm(B) does not bias β_AB. It only makes γ uninterpretable, and γ is a nuisance parameter. **No change needed** there. The fixes below concern precision of specification, the estimand, calibration and benchmark inference.

---

1. **[must] Specify df, the design and missing values exactly.** Change §5.1 step 3 to df_c = n − rank(D_c) − 2, where the rank comes from pivoted QR and includes the intercept. df_c differs by chromosome, because acrocentric chromosomes drop one arm. All pairs must share one complete-case design: give any S1 line missing MSIScore an indicator column, never drop it per pair. Targets with masked CN (311 genes) cannot be cis-adjusted. Exclude them from the universe and state this. Add a unit test: for 1,000 random pairs, t from the GEMM path equals the statsmodels OLS t to within 1e-6, including the x_B term.

2. **[should] Use a cis basis for B instead of a single linear x_B.** Homozygous deletion of B sets y_B ≈ 0, and a linear γ·x_B fits that badly. Use W_B = [x_B, 1(relCN_B < 0.25)], with an optional amplification hinge. This is still closed form. Orthonormalise each target's residualised W̃_B into Q_B (n×k). Then
⟨x̃_A, ỹ_B⟩_{·W} = x̃_A'ỹ_B − (Q_B'x̃_A)'(Q_B'ỹ_B), and ‖x̃_A‖²_{·W} = ‖x̃_A‖² − ‖Q_B'x̃_A‖².
The terms Q_B'x̃_A take k GEMMs (G_A×G_B each), so the cost is k+1 GEMMs instead of 2. Set df = n − rank(D_c) − 1 − k.

3. **[must] Pre-list arm-handling variants and state the estimand.** LOCO is the right default against cross-chromosome arm confounding, and including arm(B) is correct: it lies on the path arm(A)~arm(B)→y_B. But excluding all of chromosome c means M1 gives every gene on c the trans signature of its arm, and arm events dominate CN variance (§7 of the data profile). §1 should say that M1's β_AB is the effect of "A's dosage, including the arm-level dosage it shares". §8 allows tuning only among pre-listed variants, so add these to §5.4 now:
   - (a) **Leave-one-arm-out:** drop only arm(A) and adjust for the opposite arm of c. This separates arm-specific from whole-chromosome signal.
   - (b) **Focal:** also adjust for arm(A), so x̃_A becomes the focal deviation. Report the variance retained, var(x̃_A | +arm(A)) / var(x̃_A). In this variant, drop regulators where it is < 0.05.
   - (c) **LOCO CN-PCs:** the top 10–20 PCs of arm-adjusted CN on chromosomes ≠ c, to capture sub-arm co-selection. CRC 20q13 is an example.
   
   Each variant changes only D_c, so the per-chromosome (or per-arm) QR structure is unchanged; the arm variants need at most 41 groups instead of 22. **[could]** PureCN CIN, aneuploidy and LoH include chromosome c. A LOCO aneuploidy score (the number of derived arms ∉ c with |CN| > 0.3) would be cleaner.

4. **[should] Pre-list a small grid for CN transforms and outliers, and add influence diagnostics.** The floor 0.05 puts homozygous deletions at log2 = −4.3. Typical one-copy losses sit near −0.6 and the per-gene SD is about 0.36. With 0.5/99.5% winsorisation at n ≈ 1,037, only about 5 lines per tail are clipped, so 6–20 homdel or amplicon lines keep enormous leverage. These are informative natural knockouts, so cap their leverage rather than delete them. Grid: floor ∈ {0.05, 0.25} × winsor ∈ {0.5%, 1%}, plus INT of both x and y. INT of y matters for zero-inflated low-expression targets. For each top-K pair, report the maximum single-line share of Σx̃ỹ and the sign stability when one lineage is left out. These diagnostics are label-free, so they may also be used as filters.

5. **[must] Calibrate before BH, and define the FDR family.** The Gaussian t tails are not trustworthy at p ~ 1e-9 when both x and y are heavy-tailed, and residual structure inflates the bulk. Some context: Bonferroni at m = 1.4×10⁸ needs |z| ≈ 6.3, which is |r| ≈ 0.21 at df ≈ 935. A trans r is roughly r_cis × r_trans, so expect few honest calls.
   - (i) **Fix the NC1 permutation.** Permute the rows of the *residualised regulator matrix* X̃_A within lineage strata, jointly across all regulator genes. Keep x̃_B and Ỹ fixed. This is the residual-permutation scheme used in tensorQTL. Permuting raw CN instead decouples CN from the CN-derived covariates (CIN, ploidy, arms) and also destroys the cis adjustment for B.
   - (ii) Report λ_GC (median χ² / 0.455) and QQ plots per chromosome, for real and NC1 scores.
   - (iii) Use an empirical null. Fit Efron central matching per regulator chromosome, z′ = (z − μ̂₀)/σ̂₀, and apply BH to z′. Report BY as the conservative bound. Do not use a per-regulator null: genuine hub regulators affect up to about 48% of features and would break central matching.
   - (iv) The family is all analysed regulators × measured trans targets. Also report how many distinct CN segments are called, which is the effective number of regulators.

6. **[must] Correct the IV inference and the sign rule (§5.2).** The reduced-form t is the Anderson–Rubin test of β_IV = 0. It is not "identical" to the 2SLS Wald t, which is unreliable when instruments are weak. Report β_IV with an AR/Fieller interval: {b : t²(y_B − b·y_A on x_A | D_c, x_B) ≤ χ²₁,₀.₉₅}, which reduces to a quadratic in b. Gate on first-stage F = t²_cis ≥ 10, with τ_cis tuned at or above that. Co-segment genes violate the exclusion restriction, so label β_IV as "per unit A mRNA *if A is the sole mediator*". The knockdown-predicted sign is −sign(β_AB · β_cis,A). If anti-scaling regulators (about 5%) survive gating, using sign(β_AB) alone gets the concordance wrong for them.

7. **[must] Rename and rebuild the "mediation evidence" score (§5.3a).** The correlation of y_B with the non-CN residual of y_A is **co-expression, not mediation**. It uses exactly the variation the instrument exists to avoid, so it is confounded by shared programmes, cell state, common upstream regulators and reverse causation B→A. "Sum of z when signs agree" is discontinuous and has no stated null. Change it as follows:
   - Call it z_co ("co-expression support").
   - Combine by weighted Stouffer: z_comb = (z_RF + w·z_co)/√(1+w²), with w ∈ {0, 0.25, 0.5, 1} pre-listed.
   - Under H₀ the two statistics are independent, because the residual of y_A on x̃_A is orthogonal to x̃_A. So z_comb ~ N(0,1).
   - Accept w > 0 only if it beats both M1 and B2 on the TRAIN CMH OR.
   - Never condition z_RF on y_A to "test mediation": y_A is a collider between x_A and the confounders of y_A and y_B.

8. **[should] Make the attribution in §5.3b a single-effect posterior.** Replace the ad-hoc softmax with a SuSiE L=1 posterior built from summary statistics, over the genes g in A's window (same arm, ±10 Mb), using the same FWL z-scores:
ABF_g = √(1−ρ)·exp(ρ z²_{gB}/2), with ρ = σ₀²/(σ₀² + 1/df_c) and σ₀ ≈ 0.05 on the r scale;
α_{g|B} = π_g ABF_g / Σ_h π_h ABF_h, with prior π_g ∝ exp(λ·t_cis,g)·1(t_cis,g > 0).
Score: z_AB·α_{A|B}, or the log posterior odds. Be honest that when segment CN is nearly identical, the ABFs are equal and α ≈ π: the cis prior does all the work. Only breakpoint variation gives likelihood-based resolution. λ is tuned.

9. **[must] Make benchmark significance respect dependence (§6).** The Fisher, CMH and binomial p-values all treat pairs as independent. That fails here: the top 1% of regulators carry 55% of the hits, co-segment regulators share their score rows, and target hubs exist. Keep Fisher OR and MH OR as point estimates, but change the inference:
   - (a) **p-value from regulator-label permutation.** Reassign the Perturb-seq label rows among universe regulators within strata (own gene expressed × CN-variability tercile), and recompute the top-K OR, 1,000 times. This keeps the score structure and target hubness intact, so it tests pair-specific information. Report OR / median null OR alongside it.
   - (b) **CIs from a two-way (pigeonhole) bootstrap** (Owen 2007) that resamples regulator clusters and target clusters independently. A 5 Mb regulator bin is shorter than the arm-level CN correlation length, so use regulator arms as the primary clusters (about 39) and 5 Mb bins as a sensitivity analysis. Use target 5 Mb bins as the target clusters.
   - (c) Compare methods by a **paired** bootstrap of log-OR differences on the same resamples.
   - The MH estimator itself is fine with sparse strata (strata with no top-K pair drop out). **No change needed** there.

10. **[must] Fix the sign-concordance test.** The definition is correct: concordant iff sign(β) = −sign(LFC). Use the cis-corrected sign from item 6. The binomial test against 0.5 is wrong for two reasons. First, the base rates are not 0.5 (47.8% of HCT116 hits are up, and the sign of β is skewed by target hubs). The independence null is c₀ = q(1−u) + (1−q)u, with q = P(β>0) and u = P(LFC>0) among overlapping pairs. Second, pairs cluster by regulator. Change it to:
    - report the 2×2 sign OR and the concordance rate minus c₀;
    - test by permuting LFC signs *within regulator* among overlapping pairs, or by a regulator-cluster bootstrap;
    - also report regulator-weighted concordance, the mean of per-regulator rates over regulators with ≥ 5 overlaps;
    - compare the sign-concordant OR with a null in which prediction signs are randomly flipped within regulator.

11. **[should] Pre-declare the tuning objective as the TRAIN CMH OR at top-K,** with the pooled OR as a non-inferiority check. Cis gating and hub effects can raise the pooled OR just by ranking expressed regulators first: those average 51.0 trans hits against 5.5 for non-expressed ones. That is regulator selection, not pair knowledge. B3 shows how big this effect is, and the CMH metric removes it.

12. **[should] Restrict the GLS/LMM variant (§5.4).** Build the CN kinship LOCO, excluding chromosome c, to avoid proximal contamination, as in GWAS LMMs. Use one δ = σ²_g/σ²_e per chromosome, or a few δ bins, estimated by REML on a random sample of 500 targets. This keeps a single whitening per design and leaves df unchanged.

13. **[could] Add a heteroskedasticity-robust t as a cheap check.** For the version without x_B: HC0 Var(β̂) = Σᵢ x̃ᵢ² eᵢ² / (Σ x̃ᵢ²)². Expanding eᵢ = ỹᵢ − β̂x̃ᵢ gives Σx̃²ỹ² − 2β̂Σx̃³ỹ + β̂²Σx̃⁴. The first two sums are the GEMMs (X̃∘X̃)'(Ỹ∘Ỹ) and (X̃∘X̃∘X̃)'Ỹ, so all pairs cost two extra GEMMs. Use it for calibration diagnostics on NC1, not for ranking.

14. **[could] Deprioritise the lineage-stratified meta-analysis.** Pooled lineage fixed effects already give a within-lineage estimate. Small lineages (n = 20–40, with p ≈ 100 nuisance columns pooled) add little, and the variance gets worse. Run it last, if at all.

**No change needed:** the recursion formula and its single extra GEMM; nuisance collinearity of x_B with arm(B); ranking by |t| within a fixed universe (df_c differs by ≤ 2); per-regulator AUROC; p ≈ 100 nuisance columns against n ≈ 1,037; exclusion of expression PCs.
