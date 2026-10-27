# Round 3 review: fresh external PI (DESIGN.md v3)

Overall, the design is sound, covers the brief, and can be run once the items below are fixed. Five items are internal contradictions or undefined metrics, and each would put a wrong number into a unit test or the final run.

## Must

1. **Sign rule contradiction (§6 Sign vs §1).** §1 uses the shrunk s_A, but §6 uses raw sign(β_cis).
   *Change in §6:* "The predicted knockdown sign is −sign(β_AB)·s_A, with s_A as defined in §1. For the hinge variant, β_AB is the loss slope." *Add to §1:* "t_cis is the G1 cis t (§5.3). If A has no expression column, t_cis = 0."

2. **c₀ is undefined, and wrong under the natural reading (§6 Sign).** This was raised as round-2 ext-PI #6 but never dispositioned.
   *Change:* "q = fraction of overlapping pairs whose predicted knockdown sign is down; u = fraction with LFC < 0; c₀ = qu + (1−q)(1−u)." Pin this in a unit test.

3. **The ranked score is never defined (§5–§6).** It is unclear how t, w_G1, w_G2 and α combine, and whether ranking uses the signed or the absolute score.
   *Add §5.6:* "S_AB = t_AB · w_G1(A) · w_G2(A) · α_{A|B}. Inactive components = 1. When A3 is active, w_G2 enters only through π_g, not again as a multiplier. If Σπ = 0 in a segment, π is uniform. Direction-agnostic ranking uses |S_AB|; the predicted sign comes from §1."

4. **Resampling and permutation cannot be computed as written (§6 CIs, §8.2).** "Jackknife … with target 5 Mb bins resampled within each replicate" mixes two schemes. It is also not stated whether K is recomputed per replicate or per permuted pair set.
   *Change:* "In all jackknife, bootstrap and permutation replicates, the top-K score threshold is fixed at its full-split value. Replicate ORs are computed from per-(regulator arm × target 5 Mb bin) 2×2 counts. The primary CI is the delete-one-arm jackknife only (targets fixed). The two-way bootstrap is a sensitivity analysis. Permutations reassign labels only among regulators of the same split and cell line." Also state "K = 0.1% of the universe pairs of the evaluated split and cell line."

5. **Round-2 items silently dropped, and the TP53 release floats.** Ext-PI #6, #7 and #9 have no row in §13. §2 and §11 (T1) still say "latest programmatically available release", which contradicts the 24Q4 decision in §4 and breaks reproducibility.
   *Change:* in §2 and §11, write "DepMap 24Q4 figshare mutation matrices (pinned; checksum in DATA.md)". Add §13 rows for #6 (accepted, item 2 above), #7 and #9 (see items 8 and 13 below).

## Should

6. **Terminology clash.** PIE's "target" is our regulator A. *Add to §6:* "PIE `target` = perturbed gene = regulator A; PIE `feature` = measured gene = target B. 'BH within target' means BH within each perturbation." In the §6 universe, replace "A is a target in that line's DE table" with "A is a perturbation in that line's DE table".

7. **Adoption rule details (§8.2).** *Add:* "The 2-fold split alternates TRAIN arms within the §6 sorted order (seed 20261027) and is committed in FROZEN.sha256. If several options in a tier pass, adopt the one with the largest TRAIN Δlog-OR. The CMH guard compares against the current incumbent, not M1."

8. **Ablations relative to the selection (§8.4, round-2 ext-PI #9).** *Change:* "selected minus G2; selected minus G1; selected with the other arm handling (LOCO↔focal)." The §5.5 sensitivity analyses run inside the same single final run, listed in `selected.yaml`, and are reported only.

9. **Arm CN is not defined in the design.** *Add to §4:* "arm_k = median of x over genes on arm k (UCSC GRCh38 band at gene midpoint; 41 arms with ≥ 20 genes), computed with the active CN transform." Burden thresholds apply to this quantity.

10. **The MVP call-set fallback is undefined (§12).** *Change:* "Until Efron is built, the call set is {|t| ≥ t*}, with t* the smallest threshold where (#NC1 pairs ≥ t*/n_perm)/#real pairs ≥ t* ≤ 0.05."

11. **A3 scaling and notation (§5.4).** z_co must be on the z scale. *Change:* "z_co = atanh(r)·√(df − 1)." Rename the ABF shrinkage to "λ" (ρ is already the segment threshold) and the co-expression exponent to "w_co".

12. **Regulators with CN but no expression column** (≈ 780 CN-only genes). G1, G2 and B2 are undefined for them. *Add:* "These regulators fail G1 (t_cis = 0, s_A = +1), are 'not expressed' for G2, and B2 is set to its minimum. Their count is logged in the manifest."

13. **TP53 for chr17 regulators (round-2 #7).** Keep TP53 as a covariate: mutation is upstream of 17p LOH, so it is a confounder, not a mediator. Record this in §13 and add "TP53 covariate dropped for chr17 regulators" as a sensitivity analysis.

## Could

14. Permutation and NC1 seeds: "Seeds are derived via `numpy.random.SeedSequence(20261027).spawn(n)`; the spawn index is logged."
15. §12: "At each milestone: decision note and two-paragraph status in RESEARCH_TREE.md; prune LESSONS.md; ask a fresh agent 'what should we stop doing / what worked well'." This makes those brief duties explicit.
16. Report the TRAIN–TEST regulator CN correlation (p/q arms of the same chromosome) as a leakage diagnostic.
17. The universe is not quite "literal to the brief": Y, PAR, MT and unmapped genes are excluded. Call these operational exclusions and report their counts.

## Brief compliance

Everything is covered: cis CN with justified covariates, trans-only pairs, the universe, both overlap types (once items 1–2 are fixed), perturbed genes held out before tuning, runtime, seeds, DATA.md, and the research tree. I found no other scientific must-fix.
