# LITERATURE.md: CN-driven trans effects in DepMap, benchmarked against CRISPRi Perturb-seq

Compiled 2026-10-27 by the literature agent. I checked every reference below on the web in this session (title, authors, venue, DOI/URL). Where I could see only an abstract or a secondary summary, or where sources disagreed, the text says so. Numbers come from the cited papers. Anything I inferred is labelled *inference*.

---

## 0. Summary of the most decision-relevant numbers

| Quantity | Value | Source |
|---|---|---|
| Genes whose CN correlates positively with mRNA (tumours) | 99% of abundantly expressed genes | Fehrmann 2015 [F1] |
| mRNA change for a gained / lost arm (367 CCLE lines) | about +21% / −24% (protein +12% / −9%) | Schukken & Sheltzer 2022 [S1] |
| Genes on gained arms whose RNA "scales" | 59% (5% anti-scale) | [S1] |
| CNAs with ≥50 trans-affected genes, breast tumours | 68% (mRNA), 13% (protein) | Mertins 2016 [M1] |
| Share of a CN-regulator's trans effects matching its LINCS knockdown | about 17% on average (10 of 502 candidates significant) | [M1] |
| Median on-target CRISPRi knockdown, X-Atlas/Orion | 75.4% HCT116, 51.5% HEK293T | Huang 2025 [H1] |
| Perturbations with a discernible transcriptomic phenotype | 48.4% HCT116, 40.4% HEK293T (classifier, q<0.05) | [H1] |
| Perturbations with any downstream DE gene, K562 | 41%; 3.1% of ordered gene pairs show an effect | Aguirre 2025 [A1] |
| Typical vs essential-gene perturbation footprint | about 45 vs >500 affected genes | Nadig 2025 [N1] |
| K562 vs RPE1 similarity of perturbation profiles | median r = 0.23 (K562 day 6 vs day 8: 0.50) | Replogle 2022 [R1] |
| Best published observational-to-Perturb-seq enrichment (cis-trans eQTL colocalisation) | OR 1.9 (K562), 12.4 (HCT116); sign concordance 68% (n=231) and 71% (n=41) | Warmerdam 2026 [W1] |
| trans-eQTL links vs matched-context Perturb-seq | no enrichment above background | Yao 2023 [Y1] |

---

## 1. CNAs as natural perturbations of trans expression and protein

**Cis dosage effects are near-universal at the mRNA level but attenuated.**
- Fehrmann et al. [F1] analysed 77,840 expression profiles. Copy number correlated positively with expression for 99% of abundantly expressed genes. CN is therefore one of the dominant sources of expression covariance in cancer samples.
- Bhattacharya et al. [B1] reanalysed more than 34,000 profiles to quantify transcript-level buffering. They introduced TACNA, which uses independent component analysis to separate CNA-driven transcriptional effects from non-genetic factors.
- Schukken & Sheltzer [S1] used 367 CCLE lines with matched CN, RNA and protein. Arm calls were made relative to basal ploidy, following Cohen-Sharir [C1]. For genes on gained arms, 59% scaled at the RNA level and about 5% anti-scaled. RNA changed by about +21% on gain and −24% on loss. Only about 31% (gain) and 24% (loss) of proteins scaled. Ribosomal proteins, CORUM complex subunits, spliceosome and cell-cycle genes were the most buffered at the protein level. Extracellular and adhesion genes were enriched among anti-scaling genes. A 2026 follow-up [S2] used machine learning and named gene dependency, complex membership, haploinsufficiency and mRNA decay as predictors of buffering (abstract-level verification only).
- Gonçalves et al. [G1] used 282 CPTAC/TCGA tumours. CN effects on protein were attenuated post-transcriptionally for 23–33% of proteins, mostly complex members. They also found 48 "rate-limiting" complex interactions, where a subunit's CN drives its partners' protein but not their mRNA (e.g. AP3B1→AP3M1, GTF2E2→GTF2E1).

**Trans effects of CNAs are real and have been checked against knockdown signatures.**
- Mertins et al. [M1] (77 breast tumours) correlated each CNA gene with every mRNA and protein (Pearson, BH FDR<0.05). Trans effects were common for mRNA (68% of tested CNA genes had ≥50 trans-affected genes) and rarer for protein (13%).
  - Of 502 CNA genes queried against 3,797 LINCS knockdown signatures (connectivity score >|90|, permutation FDR<0.05), 10 were significant candidate regulators.
  - On average about 17% of a candidate's trans effects matched its knockdown signature. Examples: ERBB2, SKP1 and CETN3 on 5q (linked to EGFR), and FBXO7.
  - This is the closest precedent to our benchmark. It shows that most CN-trans associations are not recapitulated by knocking down the gene in the region.
- Li et al. (CPTAC pan-cancer, 1,064 tumours) [L1] separated cis and trans effects of drivers. Example: RICTOR gain tracks MAPKAP1 and mTOR protein.
- Goh et al. [G2] treated CN as "genomic randomisation", like Mendelian randomisation. They required three conditions:
  1. regulator expression tracks its own CN;
  2. target expression tracks regulator CN;
  3. regulator expression is only weakly correlated with the target's CN, a reverse-direction check.
  - Partial correlations beat Pearson: 16/24 vs 9/24 retrospective validation outcomes. siRNA validated 13/24 predicted pairs (54%) and 5/6 ERBB2 targets.
- Aure et al. (iPAC) [A2] residualised each gene's expression on its own CN before computing trans correlations, so that shared CN between nearby genes does not inflate them. Cis CN explained 7.5% of expression variance across all genes, 11.5% for commonly aberrant genes and 46.6% for strong in-cis genes.
- Peng et al. (remMap) [P1] used sparse multivariate regression with a group penalty to find "master predictors". They found a 17q12–q21 trans hub whose amplification affects more than 30 unlinked genes.
- Akavia et al. (CONEXIC) [A3] linked CN-altered "modulators" to expression modules.
- Wang et al. [W2] mapped CN-eQTLs in radiation-hybrid panels. 99.9% of trans effects were positive (induction). This is a reminder that generic CN-load effects can dominate trans signal.

**Aneuploidy in DepMap/CCLE.**
- Cohen-Sharir et al. [C1] called arm-level CN relative to ploidy for about 1,000 CCLE lines (997 in the supplementary table) and defined an aneuploidy score.
- The legacy "CCLE arm-level CNAs" dataset holds +1/0/−1 calls from ABSOLUTE profiles (Ghandi 2019 [G3]). DepMap staff state that ongoing aneuploidy scores are now in OmicsSignatures/OmicsGlobalSignatures [D1].

**Which genes have reliable cis CN→mRNA effects?** Most expressed genes scale at the mRNA level. Evidence across [S1, G1, B1, A2] points to:
- **Stronger (instrument-like):** metabolic genes, moderately to highly expressed genes, genes with high CN variance (focal amplicons and deletions), and in-cis genes, where CN explains up to about 47% of variance.
- **Weaker:** ribosome, complex and spliceosome subunits (buffered mainly at the protein level, partly at RNA), low-expressed genes, and extracellular genes.

**Implication for our design.**
- Estimate a per-gene cis slope and its strength (partial R² or F) for every candidate regulator A, adjusting for covariates. Use it as a weight or a filter (an "instrument strength").
- A regulator whose CN does not move its own mRNA cannot plausibly act through dosage. It may still act through protein (Gonçalves), but then the CN→trans signal is weak.
- Expect CN loss to reduce mRNA by about 20–25% on average, versus 50–75% for CRISPRi. Effect sizes will differ, but direction should be preserved.
- Include Goh's reverse check: A's expression should not be explained by B's CN.
- Protein-level complex effects (Gonçalves) are largely invisible to both mRNA readouts. Do not expect complex-partner pairs to be enriched.

---

## 2. trans-association methods and how they handle confounding

**Hidden-factor correction removes real trans signal and can create false positives.**
- PEER [S3] and PCA are standard in eQTL mapping. Zhou et al. [Z1] showed across 362 synthetic and 110 real datasets that PCA matches or beats PEER, SVA and HCP, and is far faster.
- For trans effects, factor correction is dangerous:
  - Kang et al. (ICE) [K1] showed that inter-sample expression correlation creates spurious trans "hotspots". Correcting with the global sample-correlation matrix removes them. Later work notes that this also removes genuine hotspots, because a real master regulator is itself part of the global correlation structure.
  - Fusi, Stegle & Lawrence (PANAMA) [F2] addressed this by learning hidden factors jointly with the effects of prominent regulators. They reported more trans regulators with better replication.
  - Dahl et al. [D2] showed that adjusting for PCs of molecular phenotypes induces replicating false positives, with more than 10-fold false-positive-rate inflation when covariate and confounder share about 10% variance. This is a collider-type bias.
- In our setting, Fehrmann [F1] means the top expression PCs partly *are* CN (arm-level) effects. Regressing them out removes exactly the signal we want and can induce collider bias. *Inference.*

**Gene-based and cis-mediated approaches.**
- GBAT [L2] builds leave-one-out cross-validated cis-genetic predictions of each gene's expression and tests those predictions against all trans genes. This cuts the test count by at least 100-fold. Surrogate variables are estimated *conditional on each predictor*, the "supervised SV" idea, which avoids removing the regulator's own signal.
- trans-PCO [W3] aggregates targets into gene sets (WGCNA modules, pathways) and tests a PC-based omnibus statistic, which raises power for diffuse network effects.
- Mediation:
  - Pierce et al. [P2]: of 189 trans-eQTLs that were also cis-eQTLs, 39 were significantly attenuated by the cis transcript (Sobel P<1e-5). Measurement error biases complete mediation toward partial.
  - GMAC [Y2] adaptively selects confounders per mediation test.
- In eQTLGen, trans-eQTLs replicate poorly, partly because of cell-type composition confounding [V1]. eQTLGen phase 2 [W1] inferred 47,554 directed cis→trans gene pairs by colocalising cis- and trans-eQTLs (PP4≥0.8, genes ≥5 Mb apart).

**Mendelian-randomisation and IV analogues.** A's own CN is a natural instrument for A's expression. Goh et al. [G2] used this framing. Two standard IV requirements are both threatened:
- *Instrument strength.* Weak instruments bias IV estimates (Burgess & Thompson 2011 [B2]).
- *Exclusion restriction.* A's CN segment also changes every co-segmented neighbour, any of which could be the causal gene. A's arm is also correlated with arms on other chromosomes (Section 3).
- Within-segment attribution is a fine-mapping problem. SuSiE [W4] gives credible sets of causal variables among highly correlated predictors and is a natural template.

**Linear mixed models and GLS with cell-line similarity.**
- Wainberg et al. [W5] used generalised least squares on DepMap CRISPR data to account for cell-line non-independence (lineage relatedness). They proposed a calibration check: the median p-value over all pairs should be about 0.5.
- A kinship-like matrix built from *CN profiles* (not expression) captures genome-wide CN similarity. It is the analogue of a genotype GRM: it sits upstream of expression, so it is not a collider. *Inference.* With n≈1,000, GLS whitening (one Cholesky factor applied to X and Y) costs about the same as OLS.

**Implication for our design.**
- Do **not** default to regressing out expression PCs or PEER factors.
- Covariates should be either upstream of expression (CN-derived PCs or arm-level CN, ploidy, WGD, lineage, MSI, media, sequencing batch) or estimated conditional on the regulator, GBAT/PANAMA style.
- If expression PCs are used at all, test them as an ablation and check whether known positive controls survive.
- Test both regulator CN→target expression (reduced form) and a cis-mediated version: CN→A mRNA→B, as 2SLS or as a product of a cis slope and a trans slope. Gate both on instrument strength.
- Treat co-segmented genes as a fine-mapping problem rather than assigning the signal to every gene in a segment.

---

## 3. Confounding specific to CN data

**Arm-level and whole-chromosome events.**
- Taylor et al. [T1] called arms relative to ploidy (ABSOLUTE). Their aneuploidy score is the number of altered arms (0–39). In TCGA:
  - 88% of tumours carry an arm-level or whole-chromosome event (mean score 10);
  - 8q gain and 17p/8p loss each occur in about 33–35% of tumours;
  - aneuploidy correlates with ploidy (0.55) and TP53 mutation, and correlates negatively with mutation rate in MSI/POLE tumours;
  - aneuploidy correlates positively with proliferation/E2F/G2M expression *regardless of tumour type*.
- Recurrent arm events (more than 30% of tumours) are far more frequent than the most frequent focal events (below 15%).
- Arm events co-occur in a lineage-specific way. Examples: 3p loss with 3q gain in lung squamous carcinoma [T1], and 88 synergistically co-occurring arm pairs pan-cancer (Shukla 2020 [S4]).
- Davoli et al. [D3] argue that arm patterns are selected for the *cumulative* dosage of many tumour suppressors, oncogenes and essential genes on each arm. A single arm event therefore perturbs many genes at once, and co-selected arms on different chromosomes produce trans CN correlations.

**Generic aneuploidy response.** Dürrbaum et al. [D4] studied 13 engineered aneuploid lines. They found a uniform aneuploidy response (ER, Golgi and lysosome genes up; DNA replication, transcription and ribosome genes down) that is independent of *which* chromosome is gained. CN burden thus produces trans "hub" signal that is not attributable to any specific gene.

**WGD, ploidy and CIN.**
- Zack et al. [Z2] found WGD in 37% of TCGA tumours, associated with more of every other SCNA type and with TP53 mutation.
- Quinton et al. [Q1] showed WGD+ cells have distinct dependencies (SAC, replication, proteasome, KIF18A).
- DepMap's Ploidy, WGD and CIN fields come from PureCN, and calls with ploidy >5 or poor fit are filtered [D1]. The DepMap OmicsCNGeneWGS matrix is **relative to each line's ploidy**: an unaltered gene in a tetraploid line is about 1. It has been **linear, not log2,** since 24Q2. This comes from DepMap staff on the forum [D5] rather than the formal data dictionary, so it needs checking against the release used.

**MSI.** MSI-high lines (HCT116 among them) are near-diploid, have little CN variance, and have distinct expression programs. Taylor [T1] shows the aneuploidy–mutation relationship reverses in MSI/POLE tumours. *Inference:* MSI lines contribute little to CN-trans signal and mainly act as a lineage-like covariate.

**Purity.** Purity is a dominant confounder in tumours (immune infiltration; see the Taylor leukocyte analysis [T1]). It matters less for cell lines. PureCN purity in cell lines is more a measure of fit quality than of contamination. *Inference.*

**Culture conditions.** Lagziel et al. [L3] showed that culture medium is a major confounder of metabolic-gene essentiality across DepMap lines. Media, serum-free status and growth pattern are available in Model.csv and ModelCondition.csv.

**Spatial autocorrelation.** Neighbouring genes share segments, so a trans association from A to B is shared by every gene co-segmented with A (iPAC's cis-residualisation addressed a related problem [A2]). The effective number of independent regulators is the number of distinct CN segments or breakpoints, not about 18,000 genes. This affects FDR (correlated tests) and attribution.

**Implication for our design.**
- The minimal covariate set is: lineage (OncotreeLineage, possibly subtype), ploidy, WGD, an aneuploidy/CIN score, MSI, media, and RNA-seq batch (OmicsProfiles has Stranded and SequencingPlatform).
- Add a **genome-wide CN-structure adjustment that excludes A's chromosome**: arm-level CN for all *other* arms, or leave-one-chromosome-out CN PCs (the LOCO idea from GWAS). This blocks confounding by co-selected arms on other chromosomes without removing A's own signal.
- Also adjust for B's own cis CN, and for the arm-level CN of B's arm.
- Aggregate regulators into CN segments for discovery, then attribute signal within segments using cis-expression evidence or fine-mapping.
- Expect the aneuploidy-score covariate to absorb a large generic "stress hub" component.

---

## 4. Perturb-seq resources and their use as ground truth

**Replogle et al. 2022 [R1]** (K562 genome-wide, plus K562 essential and RPE1 essential; CRISPRi dCas9-KRAB):
- 9,866 expressed genes; more than 2.5M cells; median knockdown 85.5% (K562) and 91.6% (RPE1).
- 31.1% of perturbations had a significant energy-distance phenotype (controls 1.9%). 30.5% had more than 10 DEGs.
- Phenotype strength correlated with growth defect (ρ = −0.51). 86.6% of essential-gene perturbations had a phenotype.
- 7.5% of perturbations knocked down a neighbouring gene, which justifies excluding cis pairs.
- Shared responses: integrated stress response (ISR), unfolded protein response (UPR), proteasome upregulation, cholesterol biosynthesis, NF-κB after ESCRT loss.
- Some perturbations (e.g. TTK) *induce* aneuploidy within the screen. 76% of TTK-knockdown RPE1 cells became karyotypically abnormal.
- Cross-line agreement was low: median r = 0.23 for K562 vs RPE1.

**X-Atlas/Orion (Huang et al. 2025, Xaira; bioRxiv 10.1101/2025.06.11.659105) [H1]:**
- FiCS (fix-cryopreserve) Perturb-seq; dual-sgRNA library (41,780 sgRNAs) covering all 18,903 protein-coding genes; dCas9-mScarlet-ZIM3 in HCT116 and HEK293T; about 8M cells; more than 16,000 UMIs per cell.
- Median knockdown 75.4% (HCT116) and 51.5% (HEK293T).
- Discernible phenotype for 48.4% (HCT116) and 40.4% (HEK293T) of perturbations. This was measured with L1-logistic classifiers, not DE; non-targeting false rate 4.8% and 0%.
- Strong-phenotype perturbations: 1,574 (HCT116) and 1,944 (HEK293T).
- sgRNA UMI count tracks knockdown (R ≈ 0.9), and high-sgRNA cells show dose-dependent effects.
- Batch-to-batch reproducibility was higher than K562's.

**PIE processed DE tables (arcinstitute/PIE_x_atlas_orion; Verma et al. 2026, bioRxiv 10.64898/2026.10.02.756297) [P3, P4]:**
- **Method:**
  - gpudge Mann–Whitney U test per gene, per perturbation, against that cell line's non-targeting cells, on raw counts with CP10k normalisation.
  - Genes are tested only if mean CPM ≥5 in the perturbed or control group. Untested genes are "not DE".
  - BH correction is applied per perturbation. The paper defines DE as p_adj < 0.05, with **no fold-change threshold**.
  - Columns: target, feature, target_mean, ref_mean, target_ncells, ref_ncells, log2_fold_change, p_value, Ueffect, p_adj.
- **Size:** HCT116 has 8,987 perturbed genes, 165,777 control cells and 102.2M rows. HEK293T has 4,297 perturbed genes, 218,838 control cells and 53.6M rows. Licence is CC BY-NC-SA 4.0.
- **Filtering (discrepancy to check):** The paper says X-Atlas perturbations were kept with ≥25 cells, ≥30% mean knockdown and an expressed target. It also describes a generic per-context on-target filter (target <30% of control at the perturbation level, <50% per cell, ≥30 cells). The data-profiling agent should check which filter applies to the released tables, for example via target_mean/ref_mean on the self-row.
- **Statistical caveats:**
  - The tests are per cell against very large control pools, so they are highly powered. Strong perturbations will reach significance for many tiny effects.
  - Cells are treated as independent (no pseudo-replication across guides or batches), so p-values are likely optimistic. *Inference.*
  - TRADE [N1] estimates that a typical perturbation affects about 45 genes and an essential gene more than 500.
- **Context:** PIE itself uses DepMap 25Q3 as a knowledge source, and DepMap was the most valuable source for unseen perturbations [P4]. This shows that DepMap carries information predictive of these DE tables.

**Structure of the ground truth (K562) [A1].**
- Only 41% of perturbations affect any other gene.
- 3.1% of ordered pairs have an effect.
- Out-degree (regulator) is much more heavy-tailed than in-degree (target): there are "master regulators but not master targets".

**Cell-line specificity.**
- K562 and RPE1 agree poorly [R1]. A CD4+ T cell genome-scale screen reported only partial sharing of trans effects with K562 and Jurkat (abstract-level only).
- HCT116 (ACH-000971, colorectal, MSI) is in DepMap. The DepMap HEK line is parental HEK-293 (ACH-001085, "non-cancerous"), not 293T (checked in Model.csv).

**CRISPRi vs CN loss.** CRISPRi is acute (days), strong (50–75%) and transcriptional. CN loss is chronic, typically heterozygous (about −24% mRNA [S1]) and has had time for adaptation and compensation. *Inference:* sign concordance should hold for direct dosage effects but may fail for adaptive or compensatory responses.

**Implication for our design.**
- The benchmark is a strongly hub-structured graph. Its positives concentrate on essential and strong-knockdown regulators, and on generic-response targets such as ISR, proteasome and ribosome biogenesis genes.
- Compute knockdown efficiency from the self-row and report enrichment stratified by knockdown strength and regulator out-degree.
- Consider a stricter "positive" definition as a sensitivity analysis (e.g. |log2FC| ≥ threshold, or top-k per regulator), given the very high power and the lack of a fold-change cutoff.
- Analyse HCT116 and HEK293T separately. Report HCT116 as primary because it has more perturbations and better knockdown.
- Optionally use DepMap HCT116 and HEK-293 lineage contexts as a sensitivity check, without tuning on them.

---

## 5. Benchmarking observational inference against Perturb-seq

**General GRN and forecasting benchmarks.** These are mostly negative.
- BEELINE [P5]: GRN inference from scRNA-seq had moderate AUPRC and early precision.
- CausalBench (Chevalley et al.; Communications Biology 2025 per the PMC/BioStudies record S-EPMC11897147; arXiv 2210.17283) [C2] used Replogle K562/RPE1. Methods that use interventional data did **not** outperform observational-only methods, and scalability limited performance.
- Kernfeld et al. (GGRN/PEREGGRN, Genome Biology 2025) [K2]: it is "uncommon" for forecasting methods to beat simple baselines. In the 2023 preprint, uninformed baselines matched or beat methods on MAE in every case.
- Ahlmann-Eltze, Huber & Anders [A4]: no foundation or deep model beat simple linear or additive baselines.
- Systema [V2]: "systematic variation", meaning generic perturbed-vs-control differences from selection biases and confounders, inflates common metrics. Predicting perturbation-*specific* effects is much harder.
- Aguirre et al. [A1]: baseline co-expression is essentially uncorrelated with perturbation effects across all pairs, because co-expression mostly reflects co-regulation.

**trans-eQTL and CN versus perturbation, the direct precedents.**
- **Warmerdam et al. 2026 (eQTLGen phase 2, medRxiv / PMC13119349) [W1].**
  - Method: cis-trans pairs from colocalisation (sign from allelic direction), checked against Perturb-seq pairs with BH FDR ≤ 0.05.
  - Null: degree-preserving shuffles of pairings.
  - Enrichment: OR 1.9 in K562 (P=4e-8) and OR 12.4 in HCT116 (P=2e-7). The reference list cites X-Atlas/Orion; I could not see the exact dataset mapping in the text.
  - Sign concordance: 68% of 231 (K562) and 71% of 41 (HCT116) overlapping pairs.
  - This is the best available calibration for what a strong observational trans method achieves on this benchmark. The ORs come from very few, high-confidence pairs.
- **Yao et al. 2023 (compressed Perturb-seq) [Y1].** Cis-by-trans colocalised eQTL links from LPS monocytes (N=432), a context matched to their macrophage line, were *not* enriched among Perturb-seq effects (if anything slightly depleted). They *were* enriched among co-expressed pairs.
- **Mertins 2016 [M1].** About 17% of CN-trans effects of candidate regulators matched LINCS knockdowns, and only 10 of 502 regulators were significant.
- **Morris et al. 2023 (STING-seq) [M2].** Trans networks of GWAS loci arose mainly when the cis target was a TF or miRNA.
- **Ota et al. 2025 (Nature) [O1].** Combined LoF burden effects with K562 Perturb-seq to build regulator→program→trait graphs.
- I found **no published study** that benchmarks DepMap/CCLE CN→trans-expression associations against genome-wide Perturb-seq. The closest are Mertins (LINCS) and Goh (siRNA, 24 pairs). This project appears to be novel in that respect.

**Known evaluation pitfalls.**
1. **Hub and degree bias.** Regulators with large Perturb-seq footprints (essential, strongly knocked down; Replogle ρ=−0.51 with growth [R1]) and generic-response targets dominate the positives. A method that predicts many pairs for essential regulators, or for ISR and ribosome targets, gets a high pooled OR without pair-level specificity.
   - Remedies: degree-preserving nulls [W1]; Cochran–Mantel–Haenszel ORs stratified by regulator (or by regulator and target-degree bins); per-regulator AUROC/AUPRC; a "regulator-only" baseline (score = regulator out-degree proxy) and a "target-only" baseline (score = target DE frequency).
   - FLEX [R2] showed the same problem for DepMap co-essentiality benchmarks, where a few mitochondrial complexes dominated.
2. **Universe definition.** Per the project brief, the universe is perturbed in X-Atlas, CN-variable in DepMap, target measured in both, and trans. Genes below the CPM filter are "not DE" in PIE [P4]. Untested targets therefore need to be removed from the universe, not counted as negatives. *Inference*, but this follows from the table definition.
3. **Systematic variation.** Generic perturbed-vs-control shifts (cell-cycle and stress states) [V2] mean the targets most often DE in X-Atlas resemble the expression axes most associated with aneuploidy and proliferation in DepMap [T1, D4]. This can produce spurious "concordance".
4. **Sign conventions.** CN loss should match knockdown, so positive CN→B slopes should match negative knockdown log2FC on B. Report direction-agnostic and sign-concordant results separately, with the concordance rate among overlapping pairs and a binomial test [W1].
5. **Benchmark leakage via tuning.** Hold out regulators, as in the project brief. Note that degree-structured benchmarks also leak through target-side hub identity. *Inference.*

**Typical enrichments.** The literature suggests ORs of about 1–2 for broad observational association sets (co-expression, generic trans-eQTL). Higher ORs (about 10) appear only for small, high-confidence sets backed by causal structure (cis-instrumented, colocalised) and cell-type-matched data. A realistic target here is an OR of about 1.5–3 for large prediction sets, rising for the top-ranked few hundred pairs. *Inference from [W1, Y1, M1, A1].*

**Implication for our design.**
- Make the enrichment primary, as the brief requires, but always report it alongside:
  - (i) a regulator-degree-stratified (CMH) OR;
  - (ii) degree-preserving permutation p-values;
  - (iii) naive baselines: regulator out-degree prior, target hubness prior, and plain trans Pearson correlation without covariates.
- Report enrichment as a function of top-K. Report sign concordance among overlapping pairs with a binomial test.

---

## 6. Large-scale association statistics

- **Speed.**
  - With covariates residualised out of both CN and expression (Frisch–Waugh–Lovell), all pairwise partial correlations reduce to one matrix product. Matrix eQTL [S5] and tensorQTL [T2] (GPU, more than 200× faster than CPU) are built on this idea.
  - With n≈1,000 lines, about 10⁴ regulators and about 1.8×10⁴ targets, this is about 2×10¹¹ multiply-adds, minutes on a 10-core CPU. Memory for a dense 10⁴ × 1.8×10⁴ float32 matrix is about 0.7 GB. Process in chunks and keep only the top statistics. *Inference.*
  - Adjusting each target B for **its own** cis CN is a target-specific covariate. It can still be done in closed form with the partial-correlation recursion r_{AB·c_B} = (r_{AB} − r_{A c_B} r_{B c_B}) / √((1−r²_{A c_B})(1−r²_{B c_B})), using residualised matrices. This needs one extra matrix product, so it is cheap. *Inference (standard algebra).*
- **Robustness.**
  - CN distributions are heavy-tailed (amplicons, homozygous deletions). A few lines can drive a correlation.
  - Options: winsorise or cap CN, Spearman correlation or rank-based inverse normal transform (INT) of expression, and minimum counts of lines in altered states (Schukken required ≥10 lines per category [S1]).
  - Rank-based INTs do not always control type-I error (Beasley et al. 2009 [B3]). Check the null calibration empirically.
- **FDR under dependence.**
  - BH is valid under positive regression dependence (PRDS). BY controls FDR under arbitrary dependence at about a log(m) cost in power [B4].
  - Storey's q-values [S6] estimate the null proportion π₀.
  - When there are 10⁸ correlated tests, the theoretical null is often wrong. Efron's empirical null [E1] estimates N(μ, σ²) from the bulk of the z-scores.
  - CATE / RUV-4 / LEAPP [W6] adjust multiple-testing statistics for latent confounders estimated from the test statistics themselves.
  - For an enrichment benchmark, ranking quality matters more than calibrated p-values. Calibration still matters for choosing K and for honest reporting.
- **Permutation nulls.** Permute whole CN profiles across cell lines, within lineage strata. This keeps within-genome CN correlation intact and breaks the CN–expression link. It gives a null for hit counts and for the benchmark OR. A GLS median-p ≈ 0.5 check [W5] is a cheap calibration diagnostic.
- **Shrinkage and multivariate models.**
  - remMap [P1]: ℓ1 plus group penalty to find master predictors.
  - SuSiE [W4]: fine-mapping among co-segmented regulators.
  - trans-PCO [W3]: target gene-set aggregation, which gains power for diffuse effects.
  - PANAMA [F2]: hidden factors learned jointly with regulators.
  - Low-rank (reduced-rank) regression of expression on CN is a natural middle ground between per-pair tests and full multivariate models. I have not cited a specific reference for it.

**Implication for our design.**
- Use closed-form residualise-then-correlate OLS (or GLS-whitened OLS) as the workhorse, with target-cis-CN adjustment via the partial-correlation recursion.
- Run robustness variants (winsorised, rank-based) as ablations.
- Rank by |t| (or a shrunken effect) for the benchmark. Use an empirical-null-calibrated FDR or permutation for reporting.
- Keep heavier models (SuSiE within segments, IV) for a second stage on top-ranked candidates.

---

## 7. Datasets that could help, and leakage risk

| Dataset | Use | Leakage risk vs. Perturb-seq benchmark |
|---|---|---|
| OmicsCNGeneWGS (relative-to-ploidy, linear) [D5] | regulator CN, target cis CN, arm/segment CN | none |
| Expression (log2 TPM+1, protein coding) | targets, cis-expression of A | none |
| OmicsGlobalSignatures / OmicsSignatures: MSIScore (MSISensor2), Ploidy, CIN, WGD, LoHFraction, Aneuploidy (PureCN) [D1] | covariates | none |
| Legacy arm-level CNA calls (ABSOLUTE; Cohen-Sharir 997 lines) [C1, G3] | arm covariates | none. Better to derive arm CN ourselves from OmicsCNGeneWGS for full coverage and consistency |
| Model.csv / ModelCondition.csv (lineage, subtype, MSI subtype feature, OnboardedMedia, FormulationID, SerumFreeMedia, GrowthPattern) | covariates | none |
| OmicsProfiles (Stranded, SequencingPlatform, dates) | RNA batch covariates | none |
| CRISPR gene effect (Chronos [D6]) | regulator essentiality; possible prior | **indirect.** Essentiality strongly predicts Perturb-seq out-degree [R1]. Using it to rank regulators inflates the pooled OR through hub structure, not pair-specific knowledge. Report only as a declared prior, with stratified evaluation |
| CCLE proteomics (Nusinow 2020, 375 lines [N2]); Sanger ProCan (Gonçalves 2022, 949 lines, 8,498 proteins [G4]) | protein-level mediation (A CN→A protein); attenuation filters | low, but smaller sample overlap |
| Prior networks (STRING, co-essentiality, CORUM, TF–target) | candidate priors | **moderate.** Some resources ingest perturbation or co-expression data. Verify each source's inputs and exclude any built from Replogle, X-Atlas or PIE |
| Replogle 2022, X-Atlas raw data, PIE models/embeddings, LINCS knockdown signatures | — | **high.** These are the same construct as the labels. Never use for tuning; LINCS is acceptable only as an independent secondary validation |

**Implication for our design.** All covariates the brief lists are available without leakage (purity is weakly meaningful for cell lines). DepMap CRISPR essentiality and prior networks are the main tempting shortcuts. They should be kept out of the primary method, or included only as declared, separately reported variants with degree-stratified evaluation.

---

## Design implications (ranked)

1. **Control hub and degree structure in the evaluation, not just in the model.** Report the pooled OR (primary) together with a regulator-stratified CMH OR, a degree-preserving permutation null [W1], and naive regulator-hubness and target-hubness baselines. Without these, a method can "win" by predicting essential regulators and ISR/ribosome targets [R1, A1, V2, R2].
2. **Do not regress out expression PCs/PEER by default.** CN is a major driver of expression covariance [F1]. PC adjustment removes genuine trans hubs [K1] and induces replicating false positives [D2]. Prefer upstream covariates (lineage, ploidy, WGD, aneuploidy/CIN, MSI, media, RNA batch, CN-derived PCs). Use regulator-conditional factors only if needed [L2, F2].
3. **Adjust for genome-wide CN structure leaving out A's chromosome** (LOCO arm-level CN or CN PCs), plus B's own cis CN and B's arm CN. Co-selected arms across chromosomes [T1, S4, D3] are the main trans confounder the brief warns about.
4. **Gate regulators on instrument strength.** Require a significant, positive cis CN→mRNA slope for A, and use it to weight or rank. Consider a cis-mediated (2SLS or product-of-slopes) score alongside the reduced-form CN→B score [G2, P2, B2].
5. **Treat co-segmented regulators as a fine-mapping problem.** Run discovery at segment level. Attribute to genes by cis-expression strength or SuSiE-like credible sets [W4, A2]. Otherwise every gene on an amplicon inherits the same predictions and dilutes precision.
6. **Control the generic aneuploidy/proliferation axis.** Include an aneuploidy score and CIN, and consider a proliferation score as a sensitivity covariate. Otherwise generic stress and proliferation targets will look "trans-regulated" by many CN-variable regulators [D4, T1].
7. **Define benchmark positives and universe carefully.** Drop CPM-filtered (untested) targets from the universe. Verify the knockdown filter in the released tables. Add a fold-change-thresholded positive set as a sensitivity analysis. Analyse HCT116 (primary) and HEK293T separately [P3, P4, H1].
8. **Expect modest enrichment and optimise ranking, not calibration.** Literature ORs are about 1–2 for broad observational sets and about 10 only for small, causally anchored, context-matched sets [W1, Y1, M1]. Report OR versus top-K and sign concordance with binomial tests.
9. **Engineer for closed-form speed.** Residualise, then run a matrix-product partial correlation, with the target-cis-CN recursion and GLS whitening by a CN-similarity matrix as a variant. 10⁸–10⁹ tests take minutes [S5, T2]. Calibrate with lineage-stratified CN-profile permutations and an empirical null [E1, W5].
10. **Keep leakage-prone priors out of the primary method.** Exclude DepMap CRISPR essentiality, perturbation-derived networks, Perturb-seq-trained embeddings and LINCS from tuning. If used, report them as separate labelled variants evaluated on held-out regulators.

---

## References

- [A1] Aguirre M, Spence JP, Sella G, Pritchard JK (2025). Gene regulatory network structure informs the distribution of perturbation effects. *PLoS Comput Biol* 21(9):e1013387. https://doi.org/10.1371/journal.pcbi.1013387
- [A2] Aure MR, Steinfeld I, Baumbusch LO, et al. (2013). Identifying in-trans process associated genes in breast cancer by integrated analysis of copy number and expression data. *PLoS One* 8(1):e53014. https://doi.org/10.1371/journal.pone.0053014
- [A3] Akavia UD, Litvin O, Kim J, et al. (2010). An integrated approach to uncover drivers of cancer. *Cell* 143(6):1005–1017. PMC3013278.
- [A4] Ahlmann-Eltze C, Huber W, Anders S (2025). Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nat Methods* 22:1657–1661. https://doi.org/10.1038/s41592-025-02772-6
- [B1] Bhattacharya A, Bense RD, Urzúa-Traslaviña CG, de Vries EGE, van Vugt MATM, Fehrmann RSN (2020). Transcriptional effects of copy number alterations in a large set of human cancers. *Nat Commun* 11:715. https://doi.org/10.1038/s41467-020-14605-5
- [B2] Burgess S, Thompson SG (2011). Avoiding bias from weak instruments in Mendelian randomization studies. *Int J Epidemiol* 40(3):755–764. (Citation verified; full text not read.)
- [B3] Beasley TM, Erickson S, Allison DB (2009). Rank-based inverse normal transformations are increasingly used, but are they merited? *Behav Genet* 39(5):580–595. https://doi.org/10.1007/s10519-009-9281-0
- [B4] Benjamini Y, Yekutieli D (2001). The control of the false discovery rate in multiple testing under dependency. *Ann Stat* 29(4):1165–1188. https://doi.org/10.1214/aos/1013699998
- [C1] Cohen-Sharir Y, McFarland JM, Abdusamad M, et al. (2021). Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition. *Nature* 590:486–491. https://doi.org/10.1038/s41586-020-03114-6
- [C2] Chevalley M, Roohani Y, Mehrjou A, Leskovec J, Schwab P (2025). A large-scale benchmark for network inference from single-cell perturbation data. *Communications Biology* (March 2025; PMC record S-EPMC11897147; journal DOI not verified). Preprint: arXiv:2210.17283.
- [D1] DepMap forum / release notes: 24Q2 announcement (OmicsSignatures: MSIScore, Ploidy, CIN, WGD) https://forum.depmap.org/t/announcing-the-24q2-release/3312 ; 24Q4 announcement https://forum.depmap.org/t/announcing-the-24q4-release/3564 ; aneuploidy-score threads https://forum.depmap.org/t/aneuploidy-score/3883 ; arm-level CN thread https://forum.depmap.org/t/what-happened-to-the-chromosome-arm-level-copy-number-data/3927 ; depmap_omics pipeline docs https://github.com/broadinstitute/depmap_omics/blob/master/documentation/DepMap_processing_pipeline.md
- [D2] Dahl A, Guillemot V, Mefford J, Aschard H, Zaitlen N (2019). Adjusting for principal components of molecular phenotypes induces replicating false positives. *Genetics* 211(4):1179–1189. PMC6456307.
- [D3] Davoli T, Xu AW, Mengwasser KE, et al. (2013). Cumulative haploinsufficiency and triplosensitivity drive aneuploidy patterns and shape the cancer genome. *Cell* 155(4):948–962. https://doi.org/10.1016/j.cell.2013.10.011
- [D4] Dürrbaum M, Kuznetsova AY, Passerini V, et al. (2014). Unique features of the transcriptional response to model aneuploidy in human cells. *BMC Genomics* 15:139. https://doi.org/10.1186/1471-2164-15-139
- [D5] DepMap forum: OmicsCNGeneWGS is relative to ploidy https://forum.depmap.org/t/26q1-omicscngenewgs-csv-normalization/4656 ; values linear since 24Q2 https://forum.depmap.org/t/best-practices-for-interpreting-copy-number-values-in-cn-wgs-data-25q3/4581
- [D6] Dempster JM, Boyle I, Vazquez F, et al. (2021). Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biol* 22:343. https://doi.org/10.1186/s13059-021-02540-7
- [E1] Efron B (2004). Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. *J Am Stat Assoc* 99(465):96–104. https://doi.org/10.1198/016214504000000089
- [F1] Fehrmann RSN, Karjalainen JM, Krajewska M, et al. (2015). Gene expression analysis identifies global gene dosage sensitivity in cancer. *Nat Genet* 47(2):115–125. https://doi.org/10.1038/ng.3173
- [F2] Fusi N, Stegle O, Lawrence ND (2012). Joint modelling of confounding factors and prominent genetic regulators provides increased accuracy in genetical genomics studies. *PLoS Comput Biol* 8(1):e1002330. https://doi.org/10.1371/journal.pcbi.1002330
- [G1] Gonçalves E, Fragoulis A, Garcia-Alonso L, Cramer T, Saez-Rodriguez J, Beltrao P (2017). Widespread post-transcriptional attenuation of genomic copy-number variation in cancer. *Cell Syst* 5(4):386–398. https://doi.org/10.1016/j.cels.2017.08.013
- [G2] Goh XY, Newton R, Wernisch L, Fitzgerald R (2013). Testing the utility of an integrated analysis of copy number and transcriptomics datasets for inferring gene regulatory relationships. *PLoS One* 8(5):e63780. https://doi.org/10.1371/journal.pone.0063780
- [G3] Ghandi M, Huang FW, Jané-Valbuena J, et al. (2019). Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature* 569:503–508. https://doi.org/10.1038/s41586-019-1186-3
- [G4] Gonçalves E, Poulos RC, Cai Z, et al. (2022). Pan-cancer proteomic map of 949 human cell lines. *Cancer Cell* 40(8):835–849. https://doi.org/10.1016/j.ccell.2022.06.010
- [H1] Huang AC, Hsieh T-HS, Zhu J, et al. (2025). X-Atlas/Orion: Genome-wide Perturb-seq datasets via a scalable fix-cryopreserve platform for training dose-dependent biological foundation models. *bioRxiv* 2025.06.11.659105. https://doi.org/10.1101/2025.06.11.659105 (processed data: https://doi.org/10.25452/figshare.plus.29190726)
- [K1] Kang HM, Ye C, Eskin E (2008). Accurate discovery of expression quantitative trait loci under confounding from spurious and genuine regulatory hotspots. *Genetics* 180(4):1909–1925. https://doi.org/10.1534/genetics.108.094201
- [K2] Kernfeld E, Yang Y, Weinstock JS, Battle A, Cahan P (2025). A comparison of computational methods for expression forecasting. *Genome Biol* 26:388. https://doi.org/10.1186/s13059-025-03840-y (preprint 10.1101/2023.07.28.551039; full text of the journal version not read)
- [L1] Li Y, Porta-Pardo E, Tokheim C, et al.; CPTAC (2023). Pan-cancer proteogenomics connects oncogenic drivers to functional states. *Cell* 186(18):3921–3944. https://doi.org/10.1016/j.cell.2023.07.014
- [L2] Liu X, Mefford JA, Dahl A, et al. (2020). GBAT: a gene-based association test for robust detection of trans-gene regulation. *Genome Biol* 21:211. https://doi.org/10.1186/s13059-020-02120-1
- [L3] Lagziel S, Lee WD, Shlomi T (2019). Inferring cancer dependencies on metabolic genes from large-scale genetic screens. *BMC Biol* 17. https://doi.org/10.1186/s12915-019-0654-4
- [M1] Mertins P, Mani DR, Ruggles KV, et al. (2016). Proteogenomics connects somatic mutations to signalling in breast cancer. *Nature* 534:55–62. https://doi.org/10.1038/nature18003
- [M2] Morris JA, Caragine C, Daniloski Z, et al. (2023). Discovery of target genes and pathways at GWAS loci by pooled single-cell CRISPR screens. *Science* 380:eadh7699. https://doi.org/10.1126/science.adh7699
- [N1] Nadig A, Replogle JM, Pogson AN, et al. (2025). Transcriptome-wide analysis of differential expression in perturbation atlases. *Nat Genet*. https://doi.org/10.1038/s41588-025-02169-3 (figures of 45 and >500 genes taken from the abstract)
- [N2] Nusinow DP, Szpyt J, Ghandi M, et al. (2020). Quantitative proteomics of the Cancer Cell Line Encyclopedia. *Cell* 180(2):387–402. https://doi.org/10.1016/j.cell.2019.12.023
- [O1] Ota M, Spence JP, Zeng T, Dann E, Milind N, Marson A, Pritchard JK (2025). Causal modelling of gene effects from regulators to programs to traits. *Nature*. https://doi.org/10.1038/s41586-025-09866-3 (preprint 10.1101/2025.01.22.634424; volume and pages not verified)
- [P1] Peng J, Zhu J, Bergamaschi A, et al. (2010). Regularized multivariate regression for identifying master predictors with application to integrative genomics study of breast cancer. *Ann Appl Stat* 4(1):53–77. arXiv:0812.3671.
- [P2] Pierce BL, Tong L, Chen LS, et al. (2014). Mediation analysis demonstrates that trans-eQTLs are often explained by cis-mediation: a genome-wide analysis among 1,800 South Asians. *PLoS Genet* 10(12):e1004818. PMC4256471.
- [P3] Arc Institute (2026). PIE_x_atlas_orion dataset card. https://huggingface.co/datasets/arcinstitute/PIE_x_atlas_orion (CC BY-NC-SA 4.0)
- [P4] Verma R, Adduri A, Bevilacqua B, Eraslan B, Burke DP, Goodarzi H, Roohani YH (2026). PIE: Generalizing perturbation effects across unseen perturbations, contexts and datasets. *bioRxiv* 2026.10.02.756297. https://doi.org/10.64898/2026.10.02.756297 ; code https://github.com/ArcInstitute/pie
- [P5] Pratapa A, Jalihal AP, Law JN, Bharadwaj A, Murali TM (2020). Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data. *Nat Methods* 17:147–154. https://doi.org/10.1038/s41592-019-0690-6
- [Q1] Quinton RJ, DiDomizio A, Vittoria MA, et al. (2021). Whole-genome doubling confers unique genetic vulnerabilities on tumour cells. *Nature* 590:492–497. PMC7889737.
- [R1] Replogle JM, Saunders RA, Pogson AN, et al. (2022). Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell* 185(14):2559–2575. https://doi.org/10.1016/j.cell.2022.05.013
- [R2] Rahman M, Billmann M, Costanzo M, et al. (2021). A method for benchmarking genetic screens reveals a predominant mitochondrial bias. *Mol Syst Biol* 17:e10013. https://doi.org/10.15252/msb.202010013
- [S1] Schukken KM, Sheltzer JM (2022). Extensive protein dosage compensation in aneuploid human cancers. *Genome Res* 32(7):1254–1270. https://genome.cshlp.org/content/32/7/1254 (numbers taken from bioRxiv 10.1101/2021.06.18.449005)
- [S2] Heller E, Barthel M, Räschle M, Schukken K, Sheltzer J, Storchová Z (2026). Protein buffering of aneuploidy is driven by coordinated factors identified through machine learning. *Mol Syst Biol* 22:563–598. https://doi.org/10.1038/s44320-026-00187-9 (abstract only)
- [S3] Stegle O, Parts L, Piipari M, Winn J, Durbin R (2012). Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses. *Nat Protoc* 7(3):500–507. https://doi.org/10.1038/nprot.2011.457
- [S4] Shukla A, Nguyen THM, Moka SB, et al. (2020). Chromosome arm aneuploidies shape tumour evolution and drug response. *Nat Commun* 11:449. https://doi.org/10.1038/s41467-020-14286-0
- [S5] Shabalin AA (2012). Matrix eQTL: ultra fast eQTL analysis via large matrix operations. *Bioinformatics* 28(10):1353–1358. https://doi.org/10.1093/bioinformatics/bts163
- [S6] Storey JD, Tibshirani R (2003). Statistical significance for genomewide studies. *PNAS* 100(16):9440–9445. https://doi.org/10.1073/pnas.1530509100
- [T1] Taylor AM, Shih J, Ha G, et al. (2018). Genomic and functional approaches to understanding cancer aneuploidy. *Cancer Cell* 33(4):676–689. https://doi.org/10.1016/j.ccell.2018.03.007
- [T2] Taylor-Weiner A, Aguet F, Haradhvala NJ, et al. (2019). Scaling computational genomics to millions of individuals with GPUs. *Genome Biol* 20:228. https://doi.org/10.1186/s13059-019-1836-7
- [V1] Võsa U, Claringbould A, Westra H-J, et al. (2021). Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression. *Nat Genet* 53(9):1300–1310. https://doi.org/10.1038/s41588-021-00913-z
- [V2] Viñas Torné R, Wiatrak M, Piran Z, et al. (2025). Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation. *Nat Biotechnol*. https://doi.org/10.1038/s41587-025-02777-8
- [W1] Warmerdam CAR, Westra H-J, van der Graaf A, …, Franke LH, Võsa U, et al.; eQTLGen and sc-eQTLGen Consortia (2026). Trans-eQTLs reveal the architecture of human gene regulatory networks. *medRxiv* 2026.02.04.26343575. https://doi.org/10.64898/2026.02.04.26343575 ; PMC13119349 ; code https://github.com/eQTLGen/WarmerdamEtAl2026
- [W2] Wang RT, Ahn S, Park CC, Khan AH, Lange K, Smith DJ (2011). Effects of genome-wide copy number variation on expression in mammalian cells. *BMC Genomics* 12:562. https://doi.org/10.1186/1471-2164-12-562
- [W3] Wang L, Babushkin N, Liu Z, Liu X (2024). Trans-eQTL mapping in gene sets identifies network effects of genetic variants. *Cell Genomics* 4:100538. PMC11019359.
- [W4] Wang G, Sarkar A, Carbonetto P, Stephens M (2020). A simple new approach to variable selection in regression, with application to genetic fine mapping. *J R Stat Soc B* 82(5):1273–1300. https://doi.org/10.1111/rssb.12388
- [W5] Wainberg M, Kamber RA, Balsubramani A, et al. (2021). A genome-wide atlas of co-essential modules assigns function to uncharacterized genes. *Nat Genet* 53(5):638–649. https://doi.org/10.1038/s41588-021-00840-z
- [W6] Wang J, Zhao Q, Hastie T, Owen AB (2017). Confounder adjustment in multiple hypothesis testing. *Ann Stat* 45(5):1863–1894. https://doi.org/10.1214/16-AOS1511
- [Y1] Yao D, Binan L, Bezney J, et al. (2023). Scalable genetic screening for regulatory circuits using compressed Perturb-seq. *Nat Biotechnol*. https://doi.org/10.1038/s41587-023-01964-9 (trans-eQTL comparison read in the bioRxiv version 10.1101/2023.01.23.525200)
- [Y2] Yang F, Wang J, GTEx Consortium, Pierce BL, Chen LS (2017). Identifying cis-mediators for trans-eQTLs across many human tissues using genomic mediation analysis. *Genome Res* 27(11):1859–1871.
- [Z1] Zhou HJ, Li L, Li Y, Li W, Li JJ (2022). PCA outperforms popular hidden variable inference methods for molecular QTL mapping. *Genome Biol* 23:210. https://doi.org/10.1186/s13059-022-02761-4
- [Z2] Zack TI, Schumacher SE, Carter SL, et al. (2013). Pan-cancer patterns of somatic copy number alteration. *Nat Genet* 45(10):1134–1140. https://doi.org/10.1038/ng.2760

**Verification caveats.**
- Volume and page numbers for Ota 2025 and the journal DOI for CausalBench 2025 were not confirmed.
- The Kernfeld 2025 and Nadig 2025 journal texts were not read beyond the abstract and indexing records.
- The X-Atlas knockdown filter in the PIE tables has two descriptions in the PIE paper; the data agent must check it.
- Some author lists are abbreviated with "et al." where I verified only the first authors.
- Author initials for a few multi-author papers (e.g. [A3], [M2], [Q1], [S4], [Y1]) are given as indexed by PubMed/PMC and were not individually cross-checked.
