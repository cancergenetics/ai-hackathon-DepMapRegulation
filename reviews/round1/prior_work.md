# Round 1 review: prior work (LITERATURE.md and use of prior work in DESIGN.md v1)

Reviewer: prior-work scientist. Method: web checks of 14 claims (full text where reachable, abstracts or secondary sources otherwise).

## (a) Verification table

| # | Claim in LITERATURE.md | Verdict | Source |
|---|---|---|---|
| 1 | Warmerdam 2026 (eQTLGen ph.2): OR 1.9 (K562, P=4e-8), 12.4 (HCT116, P=2e-7); sign concordance 68% (n=231), 71% (n=41); degree-preserving null | **Verified from full text.** Mapping resolved: K562 = Replogle/Ota, HCT116 = X-Atlas/Orion. Extra detail: the OR compares 10,406 colocalised cis-trans positives against 10,406 degree-matched shuffled negatives, not a CN-variable-regulator universe. The HCT116 OR rests on 41 overlapping pairs, and its sign test is only P=9e-3. | https://pmc.ncbi.nlm.nih.gov/articles/PMC13119349/ |
| 2 | Huang 2025 X-Atlas/Orion: median knockdown 75.4% / 51.5%; discernible phenotype 48.4% / 40.4%; 41,780 sgRNAs, 18,903 genes, ~8M cells | **Verified.** The same paper gives 54.3% for the K562 benchmark, which is useful context. | https://www.biorxiv.org/content/10.1101/2025.06.11.659105v1.full |
| 3 | PIE dataset card: HCT116 8,987 perturbations / 102.2M rows; HEK293T 4,297 / 53.6M; Mann-Whitney; CP10k; ≥5 CPM; CC BY-NC-SA | **Verified.** The preprocessed arrays fill untested pairs with FDR=1 (placeholders); the `de/` parquet does not. | https://huggingface.co/datasets/arcinstitute/PIE_x_atlas_orion |
| 4 | PIE paper (Verma 2026) says DepMap was the most valuable knowledge source | **Unverified.** I could not find the preprint. The PIE README confirms that DepMap CRISPRGeneEffect is an input knowledge source. | https://github.com/ArcInstitute/pie |
| 5 | Aguirre 2025: 41% of perturbations affect another gene; 3.1% of ordered pairs; heavier-tailed out-degree than in-degree; co-expression uncorrelated with perturbation effects | **Verified** (5,247 perturbations; 865,719 pairs). Co-expression is uncorrelated across all pairs but correlated among pairs sharing an edge. | https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013387 |
| 6 | Fehrmann 2015: 77,840 profiles; CN correlates positively with expression for 99% of abundantly expressed genes | **Verified** (via the Broad record). The "abundantly expressed" qualifier is kept in the text. | https://www.broadinstitute.org/publications/broad6383 |
| 7 | Mertins 2016: trans effects for 68% of CNAs (mRNA) and 13% (protein); 10 of 502 candidates significant against LINCS; ~17% match | 68% and 13% **verified** (secondary source). The ≥50-gene threshold, "10/502" and "17%" are **not confirmed**. Minor error: the text lists "ERBB2, SKP1, CETN3 on 5q", but ERBB2 is on 17q (only CETN3 and SKP1 are 5q). Also 105 tumours were profiled and 77 had high-quality data. | https://tcr.amegroups.org/article/view/10135/html |
| 8 | Schukken & Sheltzer 2022: 367 lines; 59% of RNAs scale on gain, 5% anti-scale; RNA +21% / −24%; protein 31% / 24% scale | 367 lines and the qualitative result **verified**. The percentages **could not be verified** (Genome Research paywalled, bioRxiv returned 403). Treat as unverified. Heller 2026 (MSB) exists. | https://genome.cshlp.org/content/32/7/1254 |
| 9 | Kang 2008: inter-sample correlation produces spurious hotspots | **Verified.** The later claim that the correction removes real hotspots comes from follow-up papers (Joo 2014), not Kang. | https://pmc.ncbi.nlm.nih.gov/articles/PMC2600931 |
| 10 | Dahl 2019: PC adjustment creates replicating false positives; >10-fold FPR inflation at ρ²≈10% | **Verified.** | https://pmc.ncbi.nlm.nih.gov/articles/PMC6456307 |
| 11 | Yao 2023: eQTL cis-trans links not enriched among Perturb-seq effects, but enriched among co-expressed pairs | **Unverified.** A Harvard thesis chapter from the same project says concordance with eQTLs is low, which is consistent. I did not find the numbers. | https://www.biorxiv.org/content/10.1101/2023.01.23.525200.full.pdf |
| 12 | Nadig 2025 (TRADE): typical perturbation ~45 genes, essential >500 | **Unverified.** The paper exists and TRADE/TWI is described correctly. The numbers were not found. | https://www.biorxiv.org/content/10.1101/2024.07.03.601903.full.pdf |
| 13 | CausalBench: interventional methods did not beat observational ones | **Partly outdated.** This was the 2022 preprint conclusion. The challenge report says winning solutions significantly improved on baselines, and the leading entries used interventional information. Do not cite it as a settled negative. | https://proceedings.mlr.press/v275/chevalley25a.html |
| 14 | Systema (Viñas Torné 2025, Nat Biotechnol): systematic variation inflates metrics | **Verified.** | https://brbiclab.epfl.ch/projects/systema |

No hallucinated references were found. Errors are minor: the ERBB2 chromosome slip, the unverified numerics (#4, #7, #8, #11, #12) and the stale CausalBench summary. LITERATURE's "no published DepMap-CN-vs-Perturb-seq benchmark" statement is consistent with my searches, which also found no CCLE/DepMap trans CN-QTL paper from 2024-25 (absence of evidence only).

## (b) Design changes

**1. MUST: make instrument-gated, sign-consistent scoring the primary method and ablate it, not only a tuned option.**
- Rationale: the only strong enrichment in the literature (Warmerdam OR 12.4) comes from cis-instrumented pairs. In PIE HCT116, 57% of perturbed regulators are not expressed and yield ~5 hits each (DATA_PROFILE §11), so they act as a null.
- Change: the primary score is z_trans(A→B) combined with a monotone function of the cis t of A (a TWAS or product-of-slopes form, with the cis t of A as the instrument). The sign-concordant OR is computed only where cis β > 0. Report ungated M1 as the ablation.
- Warmerdam concordance was only 68-71%, so do not expect a concordance rate above ~0.7.

**2. MUST: evaluation must not treat correlated pairs as independent.**
- Segment-sharing means hundreds of co-segmented regulators carry near-identical scores against the same B, so Fisher p-values are anti-conservative.
- Report primary p-values from block permutation or bootstrap over regulator bins (or arms) and over targets, not from Fisher.
- Add a "segment-collapsed" metric: one representative regulator per (segment, B).
- Keep Fisher OR as the point estimate only.

**3. MUST: stratify for target detectability.**
- Perturb-seq power depends on target expression (≥5 CPM filter), and 55% of hits come from the top 1% of regulators.
- Add a target-mean-expression baseline (B4 is currently variance only).
- Make the CMH stratification regulator × target-expression-decile, alongside the regulator stratification.
- Without this, the OR can simply reflect expressed targets.

**4. SHOULD: split by chromosome arm (preferably chromosome), not 5 Mb bin.**
- Arm-level events give CN correlation across bins on the same arm, so train and test regulators share DepMap signal and target hubs.
- Holding out ~6-7 chromosomes costs power but removes the leakage. Alternatively use arm-blocks and show that results are insensitive.
- The split must still be frozen before any model is run.

**5. SHOULD: add a declared DepMap-only prior variant.**
- Score = w_reg(cis F, CN variance, regulator expressed) × w_tgt(target mean expression) × z_trans.
- This is legitimate, since it uses only DepMap, and is likely the largest OR gain. It is hub exploitation, so report it only next to B3/B4 and the CMH OR to show the pair-specific gain.
- Chronos essentiality can be one further labelled variant (essentiality predicts Perturb-seq footprint). Keep it out of the primary method, as DESIGN §8 already does.

**6. SHOULD: demote M3(a) residual-expression mediation.**
- Co-expression is uncorrelated with perturbation effects across pairs (Aguirre), so a conditional co-expression term risks adding noise.
- Co-segment attribution should rest primarily on cis strength and on segment-level collapse with within-segment ranking (M3b).
- Keep (a) as a tuned weight, which may be 0.

**7. SHOULD: add a regulator-conditional hidden-factor ablation.**
- The design bans expression PCs, which is defensible (Dahl, Kang).
- LITERATURE also recommends GBAT/PANAMA-style factors estimated conditional on the regulator. They are not in DESIGN.
- Add one pre-declared ablation (expression PCs computed after residualising on CN-derived covariates and the regulator's CN, k tuned on TRAIN). It will show whether the generic stress axis is dominating.

**8. COULD: add a target-side shrinkage or module variant.**
- Perturb-seq effects concentrate on a few programmes.
- Options: a regulator-level low-rank or empirical-Bayes shrinkage of t-vectors, or trans-PCO-style gene-set aggregation, which scores B by module consistency.
- It denoises per-pair t-statistics, which are the weak link given n≈1,037.

**9. COULD: leave-HCT116-out sensitivity.**
- HCT116 is in DepMap and appears among the n≈1,037 lines.
- It is one line in a thousand, but a refit without it is cheap and settles the point.

**10. COULD: fix LITERATURE.**
- Correct the ERBB2 5q wording, and label the Schukken percentages as unverified.
- Soften CausalBench, Nadig, Yao and the PIE-DepMap claims.

**No change needed:**
- LOCO arm adjustment, B's own cis CN and B's arm CN.
- FWL closed form.
- Avoiding expression PCs by default.
- HCT116 as primary and HEK293T as secondary.
- Per-regulator AUROC and CMH.
- Permutation null NC1.
- Sign convention.
- The hub baselines.

## (c) Most promising route to high enrichment

Prior work points to a "causal-chain" signal rather than raw CN-expression correlation:
1. Restrict or weight to regulators whose own mRNA clearly follows their CN (strong, positive cis t) and that are expressed. These are the regulators whose CRISPRi can have an effect.
2. Resolve co-segment ties by cis strength.
3. Use the leave-chromosome-out adjusted reduced form, so that arm-level co-selection does not dominate.
4. Favour pairs with large |t| on highly expressed targets, because Perturb-seq power is expression-dependent.
5. Optionally add a DepMap-only regulator-hubness prior.
6. Report pair-specific gain over regulator-hub and target-hub baselines.

This mirrors what delivered OR 12.4 in the eQTL analogue, which was a small, high-confidence, cis-anchored set. Realistic expectation: OR ≈ 1.5-3 at K = 0.1%, higher at K = 0.01% (OR ≈ 3-10) with the gated, expressed-only set. Tuning on the broad primary universe will be mostly flat, so the secondary (expressed-regulator) universe should be reported prominently, and the K = 0.01% tier will be noisy.
