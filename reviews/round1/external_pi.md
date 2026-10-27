# Round 1 review: external PI (DESIGN.md v1)

Overall: this is a sound, well-scoped design. The FWL closed-form engine, LOCO arm covariates, cis-CN_B conditioning, a frozen split with a code guard, and hub baselines are all right. My main concerns are three. (i) The biggest label-free levers for the brief's objective (pooled OR) are unlisted or buried. (ii) The significance and concordance tests as written would be invalid. (iii) The variant space is too wide for a hub-dominated TRAIN set, and the rule for choosing the final method is unspecified.

## Must

1. **Pre-declare the final-method selection rule and the size of the TEST run (§8.3).**
   *Change:* add §8.7: "The final method is the single configuration with the best TRAIN primary metric, averaged over a 2-fold region-block split *within* TRAIN. A more complex variant replaces a simpler one only if its improvement exceeds the TRAIN block-bootstrap 95% CI half-width. Exactly one `--final` run evaluates the final method, all §7 baselines/controls and ≤ 3 pre-named ablations together. All results are reported."
   *Rationale:* "once per final method" leaves room for an unbounded number of final methods. A hub-dominated TRAIN set (top 1% of regulators carry 55% of hits) makes naive argmax selection noisy.

2. **Make significance robust to pair dependence (§6 Metrics).**
   *Change:* the primary significance becomes a regulator-level permutation p-value: within each regulator, shuffle scores across its universe targets (≥ 1,000 permutations), which preserves regulator out-degree. Also report a degree-preserving (regulator × target-degree-bin) permutation p-value, as in W1. The Fisher p is still reported, labelled as assuming independence.
   *Rationale:* about 60k top-K pairs that share a few hundred regulators are not independent. Fisher p-values will be astronomically small and uninformative. LITERATURE implication #1 asks for degree-preserving nulls, but DESIGN omitted them.

3. **Use the correct null for sign concordance (§6).**
   *Change:* replace "binomial test vs 0.5" with a test against the expected concordance under independence, p₀ = P(β>0)·P(LFC<0) + P(β<0)·P(LFC>0), computed from the marginals of the overlapping pairs. Also report a within-regulator sign-permutation p. The sign-concordant OR keeps its current definition.
   *Rationale:* 61% of HCT116 trans positives are down-regulated, and CN–expression slopes are skewed positive by the aneuploidy axis. A 0.5 null will show "concordance" with no pair-level sign information.

4. **Add a label-free context-expression gate for regulators as a first-class, pre-listed component (§5.2).**
   *Change:* add M2b. A regulator's score is shrunk toward 0 when A is not expressed in the application context. Use DepMap HCT116 (ACH-000971) log2(TPM+1) of A, with the threshold τ_ctx tuned on TRAIN from {0.5, 1, 2}, combined with pan-cancer cis strength. Report results with and without the gate, and the CMH OR, so that regulator-level and pair-level gains can be told apart.
   *Rationale:* regulators expressed in HCT116 average 51 trans hits, against 5.5 for non-expressed ones, and 57% of universe regulators are not expressed. CN change in a silent gene has no consequence. The gate uses only DepMap data and should be the largest legitimate gain on the primary universe.

5. **Separate focal from arm-level regulator CN (§4, §5.1).**
   *Change:* add a pre-listed variant M1-focal. x_A is additionally residualised on A's own arm CN, and that arm term is kept as a separate regulator-level "arm effect" that is not attributed to individual genes. Add baseline B6: score = partial correlation of A's *arm* CN with y_B (the same for every gene on the arm).
   *Rationale:* under M1, all ~400 genes on 8q inherit near-identical arm-driven target lists. This floods the top K with non-specific pairs. B6 shows whether gene-level CN adds anything over arm dosage, which is the confounder the brief explicitly warns about.

6. **Make the split blocks CN-aware (§6 Held-out split).**
   *Change:* replace the 5 Mb bins with contiguous blocks of ≥ 20 Mb, cut at arm boundaries (centromeres). Assign one third to TEST, stratified by chromosome, and use the same blocks for the bootstrap CIs. Record the split file's hash in RESEARCH_TREE.md at the moment it is frozen.
   *Rationale:* CN residuals stay strongly correlated over tens of Mb. With 5 Mb bins, adjacent TRAIN and TEST regulators have near-identical score profiles, so tuning leaks into TEST.

## Should

7. **Rank and cap the variant list (§5.4).**
   *Change:* run the variants in this order and stop adding once a tier fails to beat the CI rule: (1) cis gate + context gate; (2) focal vs arm; (3) score normalisation (raw |t| vs per-target standardisation of t across regulators); (4) co-segment attribution with the simple softmax; (5) rank/INT; (6) S1 vs S2. Defer GLS/kinship, lineage-stratified meta-analysis, CRC context weighting and SuSiE unless tiers 1–4 are exhausted. Cap the total at ≤ 25 logged configurations.
   *Rationale:* cell lines have no kinship structure beyond lineage and patient, and lineage fixed effects plus LOCO arms already cover population structure. Every extra variant spends TRAIN signal.

8. **Add a per-target standardisation variant (§5.1).**
   *Change:* as tier 3 above. z-score each target's t across all regulators before pooled ranking.
   *Rationale:* targets on the generic proliferation/stress axis will have high |t| with thousands of regulators and dominate the pooled top K. Whether that helps or hurts against the Perturb-seq hubs is an empirical TRAIN question.

9. **Restrict mediation (M3a) to within-segment attribution (§5.3).**
   *Change:* use the non-CN expression term only to redistribute a segment's evidence across its genes, never as an additive main-effect score. The headline method must beat B2 (co-expression).
   *Rationale:* residual co-expression mostly reflects co-regulation (Aguirre 2025). Adding it as a main effect drifts away from the CN→expression estimand the brief defines.

10. **Add covariates and deduplicate samples (§3, §4).**
    *Change:* add RNA strandedness (979/745) as a covariate. Keep one model per PatientID (prefer the one with PureCN, then the default) and exclude the 14 `EngineeredModel` lines (drug-adapted, knockout, clones).
    *Rationale:* strandedness is a known batch axis and was listed in LITERATURE but dropped. Same-patient lines are pseudo-replicates that inflate the df: 58 patients contribute 128 of the matched cell lines, and 31 patients contribute 66 lines in S1.

11. **Add combined hub baseline B5 (§7).**
    *Change:* B5 = rank(B3) × rank(B4), the product of regulator and target hubness.
    *Rationale:* this is the realistic degree-only null for a pooled OR. The method must clearly beat B5 on both the pooled and CMH ORs.

12. **Keep X-chromosome regulators in the universe (§3, §6).**
    *Change:* the universe follows the brief literally, so X-linked regulators with variable CN in **female** lines are included. The tool scores them with sex-stratified CN, and Y regulators get the lowest score. The universe should not be silently shrunk.
    *Rationale:* the brief defines the universe by "variable copy number", not by autosome. Excluding X regulators is a method choice, not a universe choice.

## Could

13. Store scores as float32, not float16 (§10). float16's 3 significant digits create ties at the top-K boundary.
14. Add a loss-side slope variant (the CN effect estimated from lines below median CN). CRISPRi mimics loss, and gains may saturate or be buffered.
15. Flag cross-chromosome paralog pairs (guide or read cross-mapping) as a sensitivity exclusion.

## No change needed
Covariate rationale (lineage, ploidy/WGD, CIN/LoH, MSI, sex, media, purity assumption); excluding expression PCs; cis CN_B conditioning via the partial-correlation recursion; the trans definition; the sign convention (concordant iff sign β = −sign LFC); HCT116 as primary; positives at p_adj < 0.05 with an LFC sensitivity analysis; the runtime and compute plan; seeds; the SIM/NC/PC controls.
