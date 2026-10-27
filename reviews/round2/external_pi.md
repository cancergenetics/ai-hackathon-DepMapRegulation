# Round 2 review: external PI (DESIGN.md v2)

## A. Round-1 items

- **#1 Selection rule:** accepted. The wording is ambiguous; see item 1 below.
- **#2 Significance:** modified, and the modification is reasonable. Regulator-label permutation preserves target hubness, which within-regulator score shuffling does not. Stratifying by own-gene status means the G2 gain is not credited in the p-value; that is honest. (Could: also report the permutation without the own-status stratum, so readers see both.)
- **#3 Sign-concordance null:** accepted. However, c₀ = q(1−u) + (1−q)u leaves q and u undefined; define them (item 6).
- **#4 Context gate (G2):** accepted correctly.
- **#5 Focal vs arm, B6:** accepted correctly.
- **#6 Split:** accepted as whole-arm units. Arms are a stronger unit than the ≥ 20 Mb blocks I proposed; fine.
- **#7–#15:** accepted or reasonably modified. The dropped lineage meta-analysis is the right call.

## B. New issues in v2

### Must

1. **Make the adoption rule unambiguous and feasible (§8.2).**
   *Change:* replace "its gain exceeds the TRAIN bootstrap 95% CI half-width" with: "A variant is adopted iff (a) Δlog-OR(top-K) vs the current incumbent is > 0 in **both** TRAIN arm-folds, and (b) the one-sided 90% lower bound of the **paired** two-way bootstrap of Δlog-OR (pooled TRAIN) is > 0. The CMH guard is the same test, with Δlog-CMH-OR lower bound > −0.05."
   *Rationale:* a single method's CI half-width on a hub-dominated set will be wide enough that nothing is ever adopted. Paired differences are the right quantity, and §6 already computes them. As written, the 2-fold split is averaged, which adds no stability check.

2. **Make tuning a greedy, tier-ordered search, and move arm handling first (§5.5).**
   *Change:* "Tiers are tuned sequentially. Each tier's grid is evaluated conditional on the adopted choices of earlier tiers, with no cross-tier factorial. The order is: arm handling (4) → G1 (1) → G2 (2) → target standardisation (5) → A3 (3) → dosage shape (6) → transforms (7) → S1/hidden factors/GLS. The cap of 30 counts every evaluated configuration."
   *Rationale:* crossing the listed grids already exceeds 30 configurations: A3 alone is ρ×s×w ≥ 12. G1's cis model and A3's segments depend on the arm handling ("A's own arm handled as in the variant"), so arm handling must be fixed first, or every later tier would be retuned.

3. **Define attribution competitors independently of the benchmark (§5.4).**
   *Change:* "A segment's competitor set is all DepMap genes in the segment that pass the frozen CN-variable and expression-measured flags, whether or not they were perturbed in X-Atlas."
   *Rationale:* if posteriors are normalised only over perturbed genes, the benchmark design leaks into scores, and scores shift between cell lines for reasons unrelated to biology.

4. **Define a minimum viable pipeline and mark the deferrable items (§12).**
   *Change:* add "M4-MVP gate: M1 + G1 + G2 + B0–B6 + NC1 + SIM + the primary metric with regulator-label permutation run end-to-end on TRAIN, and a dry-run of `final_eval.py` on a *TRAIN-only mock split*, before any tier-3+ work. Deferrable, and not required for the final run: Efron z′ (fall back to BH on the theoretical t, labelled as such), AR/Fieller intervals, HC-robust t, GLS, the hidden-factor ablation, paralog/gene-class/leave-HCT116-out sensitivities. §8.4 'all §7 baselines and controls' is satisfiable by the MVP set."
   *Rationale:* v2 has accumulated about 15 analysis components. The single-shot LOCK makes a broken final run expensive. Delivery should not depend on extras.

### Should

5. **Use a stable sign rule for weak instruments (§1, §6 Sign).**
   *Change:* "sign_kd = −sign(β_AB·β_cis) only if the cis t ≤ −3; otherwise −sign(β_AB)." Add a SIM case for this.
   *Rationale:* with w_G1 = 0.5, non-gated regulators with cis t ≈ 0 would have their predicted sign flipped at random, which degrades sign concordance for no biological reason.

6. **Define c₀ (§6 Sign).**
   *Change:* "q = fraction of overlapping pairs with predicted knockdown sign = down; u = fraction with benchmark LFC < 0; c₀ = qu + (1−q)(1−u)."
   *Rationale:* the current formula, q(1−u) + (1−q)u, is correct only if q = P(β·β_cis > 0) and u = P(LFC > 0). With the natural "both are fractions down" reading it gives the *discordance* expectation. This is a sign-convention trap: define q and u explicitly and pin them in a unit test.

7. **Apply LOCO to TP53 status and make T1 non-blocking (§4, §11).**
   *Change:* drop the TP53 covariate for chr17 regulators, as is already done for burden. If DepMap mutation calls matching 26Q3 cannot be obtained programmatically, proceed without TP53 and log this in DATA.md. Pin the mutation release; never use "latest".
   *Rationale:* TP53 mutation co-occurs with 17p loss, so adjusting for it partly removes the signal of 17p regulators (including TP53 itself, a likely HCT116 hub). Floating "latest" breaks reproducibility.

8. **Compute the X variable-CN flag on sex-centred CN (§3).**
   *Change:* "For X genes, the flag uses log2 CN minus the sex-specific median before applying ≥ 20 lines with |·| ≥ 1."
   *Rationale:* I checked this. Male lines sit at a median log2 of −0.84 on X, so 716 of 756 X genes (95%) pass the current flag, against 63% of autosomal genes. The universe flag is frozen, so this needs fixing before the freeze.

9. **Specify the TEST ablations relative to the selected method (§8.4).**
   *Change:* the ablations become "selected minus G2", "selected minus G1", and "selected with the alternative arm handling (LOCO↔focal)".
   *Rationale:* the fixed list ("focal") duplicates the selected method if focal wins.

### Could

10. Arm-split p/q correlation. Whole-chromosome events make the p and q arms of a chromosome correlated across splits. Report TRAIN–TEST regulator CN correlation as a leakage diagnostic rather than changing the unit.
11. Simplify A3 if the SIM shows the per-target posteriors are prior-dominated: use a regulator-level weight π_A, which saves the per-target ABF step.

## C. Requirement coverage
All PROMPT.md requirements are covered: covariates incl. cis CN, trans-only, the literal universe, direction-agnostic and sign-concordant overlap, held-out regulators, runtime and reproducibility. There are no internal contradictions beyond items 1, 6 and 9.
