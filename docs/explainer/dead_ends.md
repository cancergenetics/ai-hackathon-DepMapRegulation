# Dead ends (appendix)

**How to read the numbers.**
- Δ is the change in log odds ratio at the top 0.1% of pairs, on HCT116 TRAIN, against the incumbent at the time.
- "Folds" gives Δ in each of the two TRAIN arm-folds.
- "lb90" is the paired 90% lower bound (jackknife).
- The adoption rule (DESIGN §8.2) needs both folds > 0 **and** lb90 > 0.

Sources: `results/m4/tier1_compare.txt`, `results/tuning/*compare*.json`, `results/tuning_log.tsv` and
RESEARCH_TREE.md. TP ids refer to [turning_points.md](turning_points.md).

## The two that taught something

1. **I1 context-detectability fusion (D6/D7, TP7).**
   - **Result:** Δ +1.57 (folds +1.37 / +1.82, lb90 +1.31), TRAIN OR 24.1. It is the largest gain of the project
     and passed the adoption rule.
   - **Why it was not selected:** the exact pair-shuffle control C1 reaches the same OR (24.5; Δ(pair) −0.02,
     CI −0.09 to 0.06). The whole gain comes from target detectability and regulator allocation. None of it comes
     from CN pair evidence.
   - **What it taught:**
     - A strict adoption rule on a pooled OR does not protect against "winning" through hubness.
     - Within-regulator controls are required.
     - The controls themselves need known-answer tests. Before the D6a loader fix, the approximate C1 reported a
       spurious pair-specific gain of +0.063 (lb90 +0.002).
   - **What happened to it:** it was kept as the separately reported S5 context layer. On TEST it scored 18.4
     against B5's 22.8 (not different).
2. **Focal arm handling (tier 1, TP6).**
   - **Result:** Δ +0.68, but the folds were −0.09 / +1.52 (lb90 −0.62), so it was rejected.
   - **On TEST:** 4.50 against 3.85, Δ +0.16 (CI −0.37 to 0.69). Still not significant.
   - **What it taught:** a large single-fold gain on a hub-dominated benchmark is usually noise. Within-regulator
     AUROC for focal t was 0.511, against 0.507 for LOCO, so removing A's own arm does not uncover gene-specific
     signal. The paired two-fold rule earned its keep here.

Near-miss worth flagging: **I3, the shrunk-IV effect score.**
- Pooled Δ was +0.18 (folds +0.10 / +0.29, lb90 −0.10), so it was rejected.
- It was the **only variant that raised the within-regulator CMH** (ΔCMH +0.42).
- It was rejected because the objective was pooled OR. It is the best candidate for a fresh hold-out (REPORT §7.6).

## Full list (one line each)

| tier / step | variant | Δlog-OR (folds; lb90) | why it failed |
|---|---|---|---|
| 1 arm handling | leave-one-arm-out | +0.04 (−0.00/+0.10; −0.07) | no gain over LOCO |
| 1 arm handling | LOCO CN-PCs (k = 20) instead of arms | −0.02 (−0.05/+0.09; −0.13) | no gain over LOCO |
| 1 arm handling | focal (own arm adjusted) | +0.68 (−0.09/+1.52; −0.62) | one-fold gain; see above |
| 3 G2 context gate | A expressed in HCT116, thr 1 or 2, w 0 or 0.5 (4 configs) | +0.05 (lb90 ≈ −0.01 to −0.04) | redundant after G1. A strong cis gate already selects expressed regulators, and hard thresholds catch only the lowest-expression bin |
| 4 A3 | ρ = 0.98 | −0.05 (−0.05/−0.03; −0.26) | segments too fine |
| 4 A3 | cis-prior scale s = 5 | +0.00 (+0.03/−0.02; −0.49) | flatter prior loses the allocation effect |
| 4 A3 | co-expression weight w_co = 0.5 | +0.28 (+0.33/+0.21; +0.02, boot lb −0.01) | failed bootstrap bound; w_co = 0.25 passed (+0.34) but lost to w_co = 0 (+0.36) |
| 5 | target standardisation (z-score each target across regulators) | −0.03 (+0.09/−0.18; −0.21) | no gain |
| 6 | hinge (loss slope only, CRISPRi mimics loss) | −0.11 (−0.25/+0.10; −0.46) | lost power; no gain |
| 7 | CN floor 2⁻⁵ | −0.01 (−0.00/−0.02; −0.07) | no effect |
| 7 | RINT of X and Y | −0.16 (−0.25/−0.01; −0.58) | worse |
| D6 post hoc | I5 LOCO compositional TPM renormalisation | +0.05 (+0.07/+0.02; −0.03) | estimand-faithful but neutral on ranking (rows r ≥ 0.99 with the originals) |
| D6 post hoc | I4 coarser segments, ρ = 0.90 | +0.21 (−0.04/+0.51; −0.32) | one-fold gain |
| D6 post hoc | I3 shrunk-IV effect score | +0.18 (+0.10/+0.29; −0.10) | lower bound < 0; only ΔCMH-positive variant (+0.42) |
| D6 post hoc | I6 lineage-half replication (conjunction t) | +0.04 (−0.01/+0.14; −0.10) | power loss cancels any specificity |
| D6 post hoc | I1 fusion γ = 1 / γ = 2 | +1.50 / +1.57 (pass) | fails C1; reported as the context layer, not selected |
| D6 post hoc | I2 graded regulator-context prior | not run | its D6 rule (must clear a within-regulator control) cannot be met by a per-regulator weight |
| design (R1) | 5 Mb-bin hold-out | — | neighbouring bins share CN (r 0.56–0.82 at 5–10 Mb); replaced by arm split (TP3) |
| design (R1) | S1 PureCN-complete sample set as primary | — | drops HCT116, 45% of MSI-H lines and 31% of samples; became a sensitivity analysis (TEST 3.40) |
| design (R1) | co-expression "mediation" score as main effect | — | it is co-expression, not mediation; confined to within-segment use, then weight 0 |
| design (R1) | lineage-stratified meta-analysis; colorectal/MSI-H context weighting | — | underdetermined (n 20–180 per lineage against p ≈ 45); CRC/MSI lines are CN-quiet |
| deferred, never run | GLS with LOCO CN kinship; Efron empirical null; AR/Fieller IV intervals | — | deferred under the minimum-viable-pipeline rule (DESIGN §12); no time spent |
| evaluation | D4 rule \|Ueffect − 0.5\| ≥ 0.1 | — | assumed a [0, 1] scale, so it kept 98.6% of positives. Ueffect is in [−1, 1]; translated to \|Ueffect\| ≥ 0.2 (47% kept) |
| evaluation | tool call set (NC1-FDR 5%) | TEST OR 1.15 | 10.7 M genome-wide calls because λ_GC ≈ 2.1 reflects real arm structure; ranking beats calling |
| evaluation | B3 regulator-hub baseline under top-K | TRAIN fold 2 −3.78 | target-agnostic, so K covers whole rows of about 4 regulators; degenerate yardstick |
| secondary line | HEK293T | TEST OR 0.49 (6 overlaps) | not a DepMap line; 83% of hits from non-expressed perturbations |
