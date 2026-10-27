# Round 2 review: statistics (DESIGN v2)

## Verification of round-1 items

These were implemented correctly, so **no change is needed**:
- the §5.1 partial-correlation formula with Q_B;
- the pivoted-QR rank and the known-answer test;
- the NC1 residual permutation;
- λ_GC, the empirical null, BH(z′) with a BY bound;
- AR/Fieller with an F gate;
- the knockdown sign −sign(β_AB·β_cis);
- the c₀ concordance null with within-regulator sign permutation;
- the LOCO burden covariate;
- Fisher p labelled as assuming independence;
- HC-robust SEs as a check only.

Defects and new problems are listed below.

---

1. **[must] Make the cis basis rank and df per target (§4, §5.1).** For most targets the indicator 1(relCN_B < 0.25) is identically zero. Orthonormalising a zero column into Q_B then divides by zero, and the fixed k = 2 overstates the lost df. When the indicator is 1 in only one or two lines it has leverage 1, which silently deletes those lines.
   - Change: use the indicator only if it is 1 in ≥ 3 lines; otherwise W_B = [x_B].
   - Set k_B = rank(W̃_B), computed per target.
   - Set df_{c,B} = n − rank(D_c) − k_B − 1.
   - Test it: one known-answer case with an all-zero indicator and one with an indicator on 3 lines.

2. **[must] Sex-adjust the variable-CN flag for X regulators (§3).** In male lines, X relCN is about 0.5, so log2 ≈ −1. With about 700 male lines, almost every X gene passes "≥ 20 lines with |log2| ≥ 1" because of sex alone, and the universe grows spuriously. Compute the flag for X genes on log2 relCN + 1(male) (or within sex, requiring ≥ 20 lines in either sex). Freeze that version. Sex in Z already handles the modelling side.

3. **[must] Correct the regulator-label permutation (§6).**
   - (a) **Universe mismatch.** When π maps A to A′ on another chromosome, the labels of A′ for targets on chr(A′) are cis for A′, and its positives on chr(A) fall outside U_A. Counting the first set as negatives deflates null positives in a chromosome-specific way. If a method scores some chromosomes' targets highly, that makes the test anti-conservative. Change: compute each permuted OR on the pairs that are trans for both A and π(A), which keeps about 91% of pairs. Do not recode them as negatives. Use p = (1 + #{OR_null ≥ OR_obs})/(1 + 1,000).
   - (b) **Interpretation.** Strata of own-gene status × CN-variability tercile do *not* preserve regulator benchmark degree. The null therefore credits regulator selection (hub regulators ranked high), not only pair knowledge. My round-1 statement "tests pair-specific information" was wrong. Reword §6 accordingly. Add a secondary, evaluation-only permutation within strata of benchmark-degree decile × own-gene status, which isolates pair-specific signal. B5 and CMH stay as guards.

4. **[must] Replace the bootstrap clusters, since the TEST split now has about 14 arms (§6).** In v2 the split unit is the arm, so TEST holds only about 14 arms (roughly one third of about 43). A regulator-arm bootstrap with 14 clusters gives intervals that are too narrow and lumpy. Change:
   - the primary CI is a delete-one-arm jackknife on log-OR with t_{G−1} critical values (G = number of TEST arms);
   - target 5 Mb bins are still resampled within each jackknife replicate (two-way);
   - the two-way bootstrap is a sensitivity analysis;
   - the paired method comparison uses the same jackknife.
   
   State in REPORT.md that TEST precision is limited by G ≈ 14.

5. **[must] Make the SuSiE L=1 attribution internally consistent (§5.4).** v2 says co-expression enters both as a "likelihood ratio" multiplying the prior and by "weighted Stouffer". Pick one. Coherent form:
   α_{g|B} ∝ π_g · ABF(z_gB) · ABF(z_co,gB)^w, with w ∈ {0, 0.25, 0.5};
   ABF(z) = √(1−ρ)·exp(ρz²/2), ρ = σ₀²/(σ₀² + 1/df), and **σ₀ = 0.05 fixed a priori, not tuned**.
   The signed score is t_AB·α_{A|B}, because the sign is needed for concordance. Singleton segments get α = 1; say so. L=1 does not remove arm-wide signal shared across segments; the focal variant does that.

6. **[should] Handle the sign when the cis effect is weak (§1, §6).** When |t_cis| < 2, sign(β_cis) is noise. A random sign then flips about half of those predictions, even though 95% of genes have positive cis. Use −sign(β_AB·β_cis) only when t_cis ≤ −2, and −sign(β_AB) otherwise. Pre-declare and unit-test this rule.

7. **[should] Re-residualise after permuting (§5.2).** Permuting X̃ within lineage keeps it orthogonal to the lineage dummies but not to ploidy, media, burden or the arms. Project the permuted X̃ onto the orthogonal complement of D_c before scoring, so NC1 uses the same df and geometry as the real fit.

8. **[should] Define the hinge variant's estimator (§5.5 item 6).** Use x⁻ = min(x_A − m_A, 0) and x⁺ = max(x_A − m_A, 0), where m_A is the median in P. The score is the t of x⁻ given x⁺, W_B and D_c, and the predicted sign comes from the x⁻ slope. Both x⁺ and W_B vary by pair, so exact FWL needs a (1+k_B+1)-dimensional Gram per pair:
   - build it from the GEMMs ⟨x^±, Q_B⟩ and ⟨x^±, ỹ_B⟩ and the per-regulator ⟨x⁻, x⁺⟩;
   - solve it with vectorised small solves;
   - cover it with the known-answer test.
   
   A "loss-only" fit must still adjust for x⁺.

9. **[should] Base the selection rule on a paired difference (§8.2).** "Gain exceeds the TRAIN bootstrap 95% CI half-width" compares a paired difference with a marginal CI and is very conservative. Instead, a variant replaces a simpler one if the lower bound of the 95% interval for Δlog-OR (paired, same resamples or jackknife) is > 0, averaged over both TRAIN folds.

10. **[should] Score NC2 before gating (§5.3, §7).** With G2 at w = 0, class-(iii) regulators (not expressed) are all set to the minimum score, so NC2 is undefined. Evaluate NC2 on the pre-gate scores of every method. More generally, report gated methods' results within each own-gene class, so the gain from regulator selection is visible.

11. **[could] Specify the regulator-conditional factor ablation (§4).** "Remove CN effects of chromosome c" is underspecified: chromosome c has more genes than lines have rank. Define it as residualising Ỹ on the top 20 PCs of chromosome-c CN, then taking the top-k PCs (k ∈ {5, 10}). Note that it is still exposed to Dahl-type collider bias.

**No change needed:** IV/AR reporting; the F-gate grid; the call-set definition; the CMH strata (MH is robust to sparse strata); G2 as label-free; the covariate dimension against n ≈ 1,450.
