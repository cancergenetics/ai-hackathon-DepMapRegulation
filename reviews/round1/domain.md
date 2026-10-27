# Round 1 review: domain (cancer genomics / CRISPRi Perturb-seq)

Reviewer remit: biological validity of DESIGN.md v1. I read PROMPT.md, DESIGN.md, LITERATURE.md and DATA_PROFILE.md (§2–11).
Overall: the estimand, the LOCO arm adjustment, the refusal to adjust for expression PCs, instrument gating and the
HCT116-primary choice are biologically sound. The main gaps are two missing upstream confounders, a mild over-adjustment
in the genome-burden covariates, and a sign rule that breaks for anti-scaling regulators. There are also several
DepMap-internal ways to align the CN signal with what an acute CRISPRi knockdown in HCT116 can actually show.

## A. Covariates

**1. [MUST] Add TP53 functional status (from DepMap somatic mutations, damaging/hotspot calls) as a covariate.**
TP53 status sits upstream of both arms of the problem. TP53-mutant lines carry far more CN alteration, WGD and CIN
(Zack 2013; Taylor 2018), and they have lower basal expression of p53 targets (CDKN1A, MDM2, ZMAT3, SESN1, RPS27L, ...).
If TP53 status is left out, every frequently altered regulator, and above all genes on 17p, picks up spurious negative
associations with the p53 programme. This is serious because the p53 response is also a dominant generic responder in
HCT116 Perturb-seq: HCT116 is TP53-WT, and knockdown of ribosomal, essential or replication genes activates p53.
Unadjusted, the pipeline would produce "concordant" hits for the wrong reason. Mutation status is a genotype, so it is
neither a mediator nor a collider. Do not use a p53-target expression score instead: it is expression-derived and would
absorb real signal. Add the mutation table to DATA.md. HEK293T is functionally p53-null (SV40 large T); note this in the
context discussion.

**2. [MUST] Add RNA-seq technical covariates: `Stranded` (979 vs 745 models) and sequencing platform (and, if available,
library-prep era).** LITERATURE.md §3 lists RNA batch in its minimal set, but DESIGN §4 leaves it out. Strandedness
changes measured TPM for genes with antisense or overlapping transcription. It also tracks onboarding era, so it
correlates with lineage mix and with CN-calling era. It is upstream of expression and costs one column.

**3. [MUST] Compute the genome-burden covariates (aneuploidy score, CIN/fraction altered, LoH fraction) leave-one-chromosome-out.**
The PureCN aneuploidy, CIN and LoH values are functions of the whole CN profile, A's chromosome included. Conditioning on
them partly conditions on x_A, so they act as a descendant of the exposure. The bias is small for a single gene but
systematic, and largest for big or frequently altered chromosomes (1q, 8q, 17, 20). LOCO arms already span the linear
part, so the only thing these covariates add is non-linear burden. Recompute that burden from the derived arm CN,
excluding chr(A), inside the existing per-chromosome QR: count of altered arms and fraction of arms with |log2| > 0.3.
Keep PureCN ploidy and WGD as given. Both are global and largely upstream, and A contributes negligibly to them.
I checked that DepMap "CIN" is PureCN CN-based, not the CIN70 expression signature, so it is not an expression-mediator
problem.

**4. [SHOULD] Document the deliberate blocking of aneuploidy-mediated paths.** If A's dosage causes CIN (CCNE1 or AURKA
gain, MAD2L1 or TP53 loss), adjusting for burden removes the A→CIN→stress-programme path. That is a real indirect effect
under the brief's total-effect estimand. It is still the right primary choice, because acute 4–7-day CRISPRi rarely
produces karyotype change (TTK is an exception, Replogle 2022). State this in §1/§4 as a justified departure from "total
effect" and keep an un-adjusted-burden variant.

**5. [SHOULD] Use finer lineage for heterogeneous OncotreeLineages.** Lymphoid, Myeloid, CNS/Brain, Esophagus/Stomach,
Soft Tissue and Bone each mix transcriptionally distinct diseases with distinct arm patterns. Use OncotreePrimaryDisease
when it has ≥ 15 lines, falling back to lineage. Add a "solid tumours only" sample-set variant: HCT116 and HEK293T are
adherent, and haematopoietic lines carry large lineage-specific programmes that lineage fixed effects only partly remove.

**6. No change needed** for: excluding expression PCs and PEER factors; keeping proliferation as sensitivity-only (it is
a mediator for MYC, CCND1, CDKN2A and RB1 dosage); sex; media; growth pattern; MSI; the purity assumption (A5); excluding
organoids; using B's cis CN plus B's arm.

## B. Mapping heterozygous CN to near-complete CRISPRi

**7. [MUST] Define the predicted knockdown sign from the IV sign, not the reduced form.** The prediction should be
sign(LFC_B | kd A) = −sign(β_AB)·sign(β_cis,A). For anti-scaling or buffered regulators (about 5% of genes;
Schukken 2022), CN loss raises A's mRNA, so sign(β_AB) alone predicts the wrong direction. With τ_cis > 0 gating the two
rules agree, but ungated variants, B1, the S2 set and the co-segment scores would all be mis-signed. Add a toy unit test
with a negative-cis regulator.

**8. [SHOULD] Add a loss-weighted, piecewise dosage variant.** CRISPRi only lowers A. DepMap's CN variance for oncogenes
comes mostly from gains and focal amplicons, and gain-side responses saturate or anti-scale more than loss-side ones
(Schukken: about −24% on loss vs +21% on gain; protein scales even less). Fit β_loss and β_gain (a hinge at log2 relCN = 0)
and score with β_loss, or a precision-weighted blend tuned on TRAIN. Homozygous deletions (CDKN2A/MTAP, SMAD4, PTEN) are
the closest DepMap analogue of a 75% knockdown. Check that winsorising at the 0.5th percentile does not clip them for
genes with ≥ 5 deep deletions. A floor at the deletion level is enough.

**9. [SHOULD] Pre-register expectations for which regulators should transfer, and report stratified results against them.**
Transfer should be best where (i) the cis instrument is strong, (ii) the gene is expressed in HCT116, and (iii) the
product's output is read out at the mRNA level: TFs and cofactors, chromatin regulators, signalling nodes with
transcriptional outputs (KRAS/MAPK, WNT, mTOR, SREBP/sterol), and RBPs or splicing and mRNA-decay factors (NMD/IR changes
are visible in 3′ counts). Transfer should be poorest for complex subunits buffered at the protein level (ribosome,
proteasome, spliceosome core), where acute strong knockdown triggers ISR/p53/proteostasis programmes that chronic
heterozygous loss is selected to avoid. CYCLOPS-type hemizygous genes (PSMC2, SF3B1, POLR2A; Nijhawan 2012) are
intermediate. HGNC `gene_group` (already in `data/external/annotation`) gives a leakage-free class annotation for this
stratification. Report OR and concordance per class on TRAIN and TEST. The point is interpretation, not tuning.

**10. [SHOULD] Add a DepMap-internal context gate: A's expression in DepMap HCT116 (ACH-000971), and in HEK-293 (ACH-001085) for the secondary line.**
A gene that is silent in HCT116 cannot have knockdown effects there. 57% of HCT116 Perturb-seq targets are
not expressed, and they average 5.5 vs 51 trans hits. Weighting the score by DepMap HCT116 expression of A (for example
TPM ≥ 1) is biology, uses no benchmark labels, and is the legitimate way to handle the non-expressed-target quirk. It is a
regulator-level prior, though, so it lifts the pooled OR through hub structure. Declare it as a separate variant, and
judge it on the CMH OR and per-regulator AUROC, which it cannot inflate. For HEK293T, HEK-293 is only a proxy.
Optionally also down-weight regulators with damaging mutations or deep deletions in HCT116 (MSI-H frameshift targets
such as TGFBR2, ACVR2A, MLH1 silencing) [COULD].

**11. [SHOULD] Replace the "colorectal / MSI-H" context weighting with heterogeneity-aware pooling.**
MSI-H and CRC lines are near-diploid and CN-quiet, so up-weighting them removes most of the instrument variance. A better
route to transferable effects is to favour effects that are consistent across contexts. Use the planned
lineage-stratified meta-analysis and down-rank pairs with high between-lineage heterogeneity (for example meta-z × (1 − I²)
or a random-effects z). A secondary variant worth testing is a TP53-WT × CN interaction, or a TP53-WT-only estimate,
because p53-dependent trans effects in HCT116 are expressed only in p53-WT backgrounds (about 30% of DepMap).

**12. [SHOULD] Add Goh's reverse check to the co-segment and mediation scores (§5.3a).** Residual (non-CN) co-expression
of A and B mostly reflects co-regulation (Aguirre 2025; Yao 2023), not A→B causation. Use mediation evidence only when
its sign agrees with the reduced form, and penalise pairs where x_B predicts y_A as strongly as x_A predicts y_B. Always
compare against B2, so that the mediation term is not just re-importing co-expression modules such as proliferation.

**13. [COULD] TPM compositional artefact.** Massive amplicons of highly expressed genes (COLO320 MYC at about 120 copies)
reduce every other gene's TPM. Winsorisation mostly handles this. A diagnostic that flags regulators whose trans
β are > 80% negative across all targets would catch residual cases.

## C. Benchmark quirks, handled without changing labels

**14. [SHOULD] Split the secondary-universe "expressed" flag into three states.** A missing self-row can mean the gene is
truly not expressed, or that it is low (below CP10k 0.05; many TFs and signalling genes) with a knockdown too complete to
test. Define (a) self-row tested and knocked down, (b) no self-row but DepMap HCT116 TPM ≥ 1, (c) neither. Use (c) as the
empirical-null set NC2. Report (b) separately. Freeze this before any model is run. It is a stratification, not a label
change.

**15. [SHOULD] Add target-degree stratification and a hub-excluded sensitivity analysis.** DESIGN stratifies CMH by
regulator only. Generic responders (ISR: ATF4, DDIT3, ASNS, TRIB3; p53 targets; histone and E2F genes; cholesterol
genes) dominate the target side, and aneuploidy and proliferation in DepMap load on the same genes (Systema bias). Add
CMH over regulator × target-degree bins, with degree computed from the benchmark for stratification only, never as a
feature. Also report the primary OR with the top 1% of regulators by benchmark out-degree removed. Neither changes the
labels. Both show whether the gains are pair-specific.

**16. [SHOULD] Give the MW-on-large-pools quirk an effect-size sensitivity based on `Ueffect`, not LFC.** With about
165k control cells, tiny shifts reach significance, and cell-state composition (fewer cycling cells after an essential
knockdown) produces broad "hits". LFC is unusable when the perturbed mean is 0. Add the sensitivity positive set
p_adj < 0.05 & |Ueffect − 0.5| ≥ δ, with δ fixed now (for example 0.05), alongside the LFC version. Keep the primary
label as released.

**17. No change needed** for: HCT116 primary and HEK293T secondary; the literal primary universe; the sign convention
itself (beyond item 7); excluding X/Y regulators; excluding cis pairs.
