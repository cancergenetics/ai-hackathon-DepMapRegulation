# Round 2 review: domain (cancer genomics / CRISPRi Perturb-seq)

## A. Round-1 dispositions

Handled correctly, no change needed: #1 TP53, #2 RNA technical, #3 LOCO burden (dropping PureCN CIN/LoH is right),
#4 blocked aneuploidy paths (now documented), #5 finer lineage, #7 IV sign, #9 class stratification, #10 G2,
#13 compositional artefact, #15 hub-excluded and target-expression CMH strata, #16 Ueffect.

The rejections are reasonable for #11 (lineage heterogeneity pooling, now that the meta-analysis is dropped) and #12
(Goh reverse check, deferred because co-expression now only redistributes weight within a segment). Two round-1 items
were only partly carried through. They are items 3 and 5 below.

## B. New or remaining issues

**1. [MUST] Shrink the cis sign in the knockdown sign rule.** sign_kd = −sign(β_AB·β_cis) uses the *estimated*
cis sign. If G1 is off or w_G1 = 0.5, weak-instrument regulators stay in the ranking. For them, sign(β̂_cis) is close
to a coin flip, while the prior is strong: 95% of genes have a positive cis correlation (DATA_PROFILE §4). Use −1 for
β_cis only when cis t ≤ −3 (pre-listed), and +1 otherwise. Without this, sign concordance is pulled toward 0.5 by
noise, not biology. The SIM negative-cis regulator should have a strong instrument, and you should add a
weak-positive-cis case to SIM.

**2. [MUST] Compute genome burden on autosomal arms only.** Relative CN of X in male lines is about 0.5 against ploidy
(log2 ≈ −1), so Xp and Xq count as "altered" in nearly every male line. Female lines carry Xi loss or gain that has
little effect on expression. The sex covariate absorbs only the mean shift. It does not absorb the interaction with
the rest of the burden, or loss of X in female lines. For X regulators, LOCO burden already excludes X, so the same
rule applies to them.

**3. [SHOULD] Refine NC2 with DepMap expression (round-1 #14, partly adopted).** Class (iii) ("own feature not
tested") includes functional regulators below 5 CPM in 3′ counting, typically TFs and signalling genes. Define NC2 as
class (iii) **and** DepMap HCT116 log2(TPM+1) < 1, and report the remainder of class (iii) separately. Freeze this
with the other flags. Evaluate NC2 on scores computed without G2, because G2 removes most of NC2 by construction, so
the null test would be empty.

**4. [SHOULD] Apply G2 inside attribution when G2 is enabled.** π_A ∝ softmax(cis t/s) tends to land on strongly
expressed passengers (for example GRB7/PGAP3 beside ERBB2). Multiply π_A by the G2 weight, so that segment signal is
not attributed to a neighbour that is silent in HCT116 while an expressed gene in the same segment is down-weighted.
With near-identical ABFs the posterior is essentially the prior, so this choice drives attribution.

**5. [SHOULD] Add a TP53-WT-only sample-set variant (round-1 #11, cheap remainder).** This is a single configuration,
like S1, with no lineage fits. p53-dependent trans effects, the major responder axis in TP53-WT HCT116, occur only in
p53-WT backgrounds (about 30% of DepMap). Report it as a sensitivity variant. HEK293T is p53-inactive (SV40 large T),
so do not expect the variant to help there. Use that as a built-in check.

**6. [SHOULD] Put the hinge at log2 relCN = 0 (copy-neutral relative to ploidy), not at the median.** On arms that
are often lost (8p, 17p, 18q) the median is below neutral, so "below the median" mixes neutral and lost lines.
Require ≥ 20 lines on each side, or else fall back to linear.

**7. [COULD] Add a regulator homozygous-deletion term in the loss-only variant.** With the CN floor at 2^-2, deep
deletions of CDKN2A, PTEN or SMAD4 have their leverage capped. Those deletions are the closest analogue to a 75%
knockdown. Adding 1(relCN_A < 0.25) as a second regulator regressor, with a 2-df score, recovers them.

**8. [COULD] Improve ploidy/WGD imputation for the about 40% without PureCN.** Imputing a binary WGD by the lineage
median is crude. Where it exists, use Sanger CMP ploidy (618 lines; r = 0.52) before the lineage median, and keep the
missingness indicator.

**9. No change needed:**
- TP53 encoding: keep it mutation-only. Do not add MDM2 amplification to the flag. That flag would be a function of
  regulator CN (chr12) and would adjust away MDM2's own trans signal.
- G2 thresholds, and HEK293 as the proxy for HEK293T.
- The residual co-expression restriction.
- The decision to block aneuploidy-mediated paths.
