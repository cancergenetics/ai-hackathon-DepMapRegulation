# Data profile (M1): DepMap 26Q3 CN/expression and the PIE X-Atlas/Orion benchmark

All numbers come from `data/processed/{depmap_profile,benchmark_profile,gene_annotation_summary,release_evidence}.json`,
which are written by `scripts/data/*.py` (run `scripts/data/run_all.sh`). Sources, versions and checksums are in `DATA.md`.
Assumptions are marked **A#**; decisions left to the PI or the design are marked **Q#**.

## 1. Release

The DepMap data are **DepMap Public 26Q3**, inferred from release-note markers; the portal md5 listing is bot-protected and could not be checked (details in DATA.md §1).
The PIE benchmark is pinned at HF commit `abfa3da0`.

## 2. DepMap samples and profile matching

| | CN (OmicsCNGeneWGS) | RNA (TPMLogp1 protein-coding) |
|---|---:|---:|
| profile rows | 1,846 | 2,314 |
| models | 1,763 | 2,130 |
| rows with `IsDefaultEntryForModel == Yes` | 1,763 (exactly 1 per model) | 2,130 (exactly 1 per model) |
| models with >1 profile row | 80 | 172 |

- **Matching (A1):** for each modality we take the single row with `IsDefaultEntryForModel == "Yes"` and join on `ModelID`.
  Every SequencingID is found in `OmicsProfiles.csv` (CN rows are all `wgs`, RNA rows all `rna`), and the default flags agree with it 100%.
  No duplicate SequencingIDs; non-default rows (other ModelConditions/re-sequencing) are ignored.
- **Matched models: 1,724** (1,509 cell lines, 215 organoids). There are 39 CN-only and 406 RNA-only models.
  In 1,688/1,724 (97.9%) matched models the CN and RNA come from the same ModelCondition; the other 36 differ (`same_mc` flag).
  Using all rows instead of defaults does not add models (still 1,724 with both).
- RNA strandedness: 979 stranded and 745 unstranded (a candidate technical covariate). The CN sequencing year ranges from 2012 to 2025.
- Sex: 922 male, 724 female, 78 unknown/missing. 14 engineered models.
- Benchmark lines in DepMap: **HCT116 = ACH-000971 (matched)**. **HEK293T is absent**; HEK293 (ACH-001085) is matched, and HEKTE (ACH-000049) has RNA only.

## 3. Genes

| | n |
|---|---:|
| CN gene columns | 19,339 (18,778 `SYMBOL (ENTREZ)` + 561 bare `ENSG…`, of which 3 are `_PAR_Y`) |
| expression gene columns | 19,218 (all `SYMBOL (ENTREZ)`, all HGNC protein-coding) |
| in both (identical column name) | **18,559** (18,560 by HGNC id) |
| CN-only | 780 = 557 without HGNC (bare Ensembl ids), 125 readthrough loci, 21 pseudogenes, ... |
| expression-only | 659 protein-coding genes absent from CN (segmental-duplication/low-mappability masking, e.g. AGAP*, ACSM2A/B, AKR1C1) |

Shared genes by chromosome class: 17,785 autosomal, 748 X, 13 Y, 13 MT.

### Gene annotation (`data/processed/gene_annotation.tsv`)

- **Mapping:**
  - All 19,437 `SYMBOL (ENTREZ)` columns map by Entrez to HGNC (complete set 2026-10-06); 11 symbols were renamed since DepMap's build (e.g. AATK→LMTK1, IL1RL1→IL33R), and the current symbol is kept in `symbol`.
  - The 561 Ensembl columns map by Ensembl id (557 have no HGNC entry).
  - No DepMap gene maps to the same HGNC id twice.
- **Coordinates:**
  - GRCh38 start/end come from Ensembl 116 for 19,967/19,998 columns. 7 more have chromosome/arm from the HGNC location only, and 24 are **unplaced**: Ensembl ids retired before Ens116, all CN-only and non-HGNC.
  - Arm and cytoband come from the UCSC band at the gene midpoint. They disagree with HGNC's location field for 0 genes (chromosome) and 0 genes (arm).
- **Chromosome classes and PAR genes:**
  - Classes: autosome / X / Y / MT / unplaced.
  - 22 pseudoautosomal (PAR) genes are flagged `is_par` and use their Ensembl X location; the 3 `_PAR_Y` CN columns are placed on Y using GRCh38 PAR offsets.
  - **A2:** trans/cis uses Ensembl chromosomes, so PAR genes count as X. Excluding PAR genes from pairs with X/Y partners is a conservative option.
- MT genes (13) have expression but **all-NaN CN**.

### Perturb-seq identifiers (`data/processed/pie_gene_map.tsv`)

`target` values are gene symbols. `feature` values are symbols, Ensembl ids (where no symbol; e.g. 3,376 in HCT116), or 10x-style `SYMBOL_ENSG…` for duplicated names.
The vocabulary has 10,007 targets and 22,182 features across both lines.

| | HCT116 targets | HEK293T targets | HCT116 features | HEK293T features |
|---|---:|---:|---:|---:|
| n | 8,987 | 4,297 | 20,348 | 18,518 |
| exact current HGNC symbol | 8,917 | 4,264 | 16,710 | 15,928 |
| Ensembl id (in Ens116 / retired) | – | – | 3,113 / 263 | 2,167 / 176 |
| rescued: HGNC previous symbol | 40 | 18 | 140 | 138 |
| rescued: Ensembl gene_name (HGNC renamed since PIE reference) | 26 | 14 | 58 | 57 |
| rescued: HGNC alias / embedded Ensembl id | 0 / 0 | 0 / 0 | 3 / 12 | 1 / 10 |
| Ensembl name without HGNC entry | 2 | 1 | 9 | 10 |
| ambiguous (multiple previous-symbol hits) | 0 | 0 | 1 (C18orf21) | 1 |
| collision (two PIE ids → one HGNC id; lower priority dropped) | 1 | 0 | 31 | 23 |
| unmapped | 1 | 0 | 8 | 7 |
| **linked to DepMap CN** | **8,766** | **4,108** | – | – |
| **linked to DepMap expression** | 8,872 | 4,195 | **14,494** | **14,052** |
| linked to both | 8,712 | 4,062 | 14,278 | 13,869 |

- **A3, rescue priority:** current HGNC symbol > unique Ensembl-116 gene_name > unique HGNC previous symbol > unique HGNC alias. Ambiguous hits are left unmapped.
- **Collisions** (38 ids in total, almost all lncRNAs) come from HGNC locus merges where an old symbol is still a separate feature in PIE's reference (e.g. LINC02910→MIR646HG).
  The higher-priority id keeps the DepMap link, so both maps stay **one-to-one** (tested).

## 4. CN value scale (verified empirically)

**`OmicsCNGeneWGS` is linear relative copy number** (not log, not absolute), normalised so that a line's
typical genome copy number is ≈ 1:
- Values are all > 0, with no negatives and no exact zeros. Homozygous deletions are encoded as ~1e-6 (e.g. CDKN2A in A549 and U87MG; Y genes in female lines).
- The genome-wide median is 1.025 (1st–99th percentile 0.45–2.09; max 260). The per-model median ranges 0.83–1.21 (median 1.03). 1.3% of values are > 2 and 0.09% are > 4.
- Comparison with 24Q4 PureCN absolute CN on 687 models where the signature is from the same profile:
  - The per-model median of absolute/relative tracks PureCN ploidy: Pearson **r = 0.976**; abs/rel ≈ 0.89 × ploidy (IQR 0.86–0.92).
  - The median per-model Pearson of relative CN vs absolute/ploidy is 0.88 (log–log 0.92).
  - So relative CN ≈ absolute CN / (≈1.1 × ploidy); ploidy itself is **not** recoverable from this file.
- Under a log2(x+1) transform the maximum would be impossible (2^260), and the 24Q2 release note confirms the linear scale.
- **Positive controls** (relative CN, percentile within 1,724 matched models):

  | gene | lines |
  |---|---|
  | ERBB2 | SKBR3 11.7, BT474 14.2, HCC1954 66.9, NCIN87 46.1, OE19 55.9 (all ≥ 98.8th pct); MCF7 0.55 |
  | MYC | COLO320 123 |
  | EGFR | A431 8.3, MDAMB468 34, HCC827 28 |
  | MDM2 | SJSA1 33 |
  | CCNE1 | HCC1569 17 |
  | CDKN2A | 28.9% of lines < 0.25 |
- **Cis CN–expression Spearman** across 18,290 shared genes: median 0.22 (IQR 0.11–0.36); 95.4% are positive and 34.8% are > 0.3. Known amplicons: CDKN2A 0.67, EGFR 0.45, MYC 0.41.
- **A4:** for variance statistics and derived arm CN we use log2(relative CN) with a floor of 2^-5. 0.09% of values fall below the floor.

## 5. Missingness

- **CN:**
  - 311 genes are NaN in *all* models (masked); every other gene is complete. So every model has exactly 1.61% missing and no model has > 5%.
  - Missingness by chromosome class: autosomal 1.3%, X 6.3%, Y 18.8%, MT 100%.
  - Y CN: in female lines the median is ~1e-6 (as expected); in male lines it is 0.29 (IQR 1e-6–0.49), i.e. frequent loss of Y.
- **Expression:** no NaN. 5 genes are zero in all models, and 7,833 genes have median log2(TPM+1) < 1 (low expression; a filter is needed in modelling).

## 6. Per-gene CN variability (1,724 matched models, log2 relative CN)

The SD of log2 CN per gene has median 0.36, IQR 0.33–0.40, minimum 0.24 (excluding all-NaN genes); the MAD median is 0.17.
**CN variation is pervasive**, so loose definitions call every gene variable:

| definition of "variable CN" | all CN genes | autosomal | HCT116 targets | HEK293T targets |
|---|---:|---:|---:|---:|
| (any non-masked) | 19,028 | 18,282 | 8,727 | 4,074 |
| SD > 0.1 / > 0.2 / ≥ 50 lines with \|log2\| > 0.3 / ≥ 10% lines with \|log2\| > 0.3 | 19,028 | 18,282 | 8,727 | 4,074 |
| SD > 0.3 | 16,807 | 16,062 | 7,713 | 3,659 |
| ≥ 5% of lines with gain > 0.3 **and** ≥ 5% with loss < −0.3 | 15,815 | 15,648 | 7,225 | 3,394 |
| ≥ 10 lines with \|log2\| > 1 (≥2-fold vs own ploidy) | 17,744 | 17,001 | 8,127 | 3,823 |
| ≥ 20 lines with \|log2\| > 1 | 13,742 | 12,999 | 6,357 | 3,041 |
| SD > 0.5 | 1,170 | 563 | 584 | 311 |

- **Q1:** the PROMPT's "variable copy number" filter has to be stricter than SD > 0.2, or it is a no-op.
  "≥ 20 lines with ≥ 2-fold change" keeps about 70% of targets. SD > 0.5 keeps 7% of targets, and 607 of its 1,170 genes are on X/Y (sex-chromosome dosage).
  Per-gene statistics are in `data/processed/depmap_cn_variability.tsv`.
- **Sex chromosomes:** relative CN of X/Y genes varies mostly with sex, so it should be excluded from regulator CN variability or modelled with sex as a covariate (**Q6**).

## 7. Arm-level CN (derived; `depmap_arm_cn_derived.parquet`)

- Arm CN is the median log2 relative CN of genes on each arm, for all 1,763 default-CN models.
- 41 arms are usable (≥ 20 genes); 21p, Yp and Yq are dropped, and acrocentric p-arms have no genes. The median is 427 genes per arm.
- 17.1% of model × arm values are gains (> 0.3) and 17.0% are losses (< −0.3).
- **Inter-chromosomal arm correlations are modest:** median |r| 0.04, 99th percentile 0.19, max 0.27, and no pair has |r| > 0.3.
  The top pairs are 18q–20q (−0.27), 18q–8p (+0.25) and 19p–21q (+0.22); these are the colorectal 18q-loss/20q-gain patterns.
  Such lineage-linked arm events are the expected trans confounders.
- Max |corr(arm CN, PureCN ploidy)| is 0.30.

## 8. Lineage, model type, media and culture

- **Lineages and subtypes:**
  - Matched models span 30 OncotreeLineages (24 with ≥ 10 models, 21 with ≥ 20), 98 primary diseases and 191 Oncotree codes; 8 models have no lineage.
  - Largest lineages: Lung 180, CNS/Brain 154, Esophagus/Stomach 141, Lymphoid 132, Bowel 102, Skin 101, Breast 98, Bone 85, Pancreas 81, Soft Tissue 79, Ovary 76, Head and Neck 68, PNS 61, Myeloid 59, Uterus 44, Kidney 39, Biliary 38, Bladder 33, Liver 25, Eye 24, Pleura 20.
  - Smaller lineages with ≥ 10 models: Cervix 19, Prostate 18, Thyroid 14.
  - Lineages with < 10 models: Ampulla of Vater 9, Testis 7, Fibroblast 4, Vulva/Vagina 2, Adrenal 2, Other 1.
- Model type: 1,509 cell lines and 215 organoids. All organoids, and only organoids, have GrowthPattern "Dome", so GrowthPattern encodes organoid status.
- **Media and culture fields** (matched models; "≥10" = levels with ≥ 10 models and the share of models in them):

  | field | levels | missing | levels ≥ 10 (share) | top |
  |---|---:|---:|---|---|
  | Model.OnboardedMedia (MF id) / Model.FormulationID (text) | 168 | 54 | 33 (79%) | RPMI+10%FBS 471, DMEM+10%FBS 166, OPAC:Ad+++ 83 |
  | ModelCondition.GrowthMedia / FormulationID (RNA condition) | 168 | 41 | 35 (81%) | RPMI+10%FBS 476, DMEM+10%FBS 163 |
  | SerumFreeMedia (Model / MC) | 2 | 54 / 41 | – | False 1,466 / True 204 |
  | Model.GrowthPattern | 6 | 0 | 5 | Adherent 1,075, Suspension 248, Dome 215, Mixed 101, Spheroid 76 |
  | ModelCondition.CellFormat | 6 | 7 | 5 | adherent 1,150, suspension 238, dome 169 |
  | PlateCoating | 3–4 | 1,642 | – | ULA 66, Laminin 12, Matrigel 4 (Model) |
  | ModelCondition.DataSource | 1 | 0 | – | BROAD (all default profiles) |

  In `ModelCondition.csv` the column `GrowthMedia` holds MF-ids and `FormulationID` holds text. They look swapped relative to their names, but they carry the same information, one-to-one.
  **Q2:** for a media covariate, collapse to the base medium (RPMI/DMEM/…) + serum-free flag, or keep the levels with ≥ 10 models.

## 9. Covariates: ploidy, WGD, MSI, CIN, LoH, purity

| covariate (source) | matched models with value | notes |
|---|---:|---|
| PureCN Ploidy, WGD, CIN, LoHFraction, Aneuploidy (24Q4) | 1,037 (60%) | cell lines only; **0 organoids**; 905 from the *same* WGS profile as the 26Q3 CN, the rest model-level (may be WES) |
| MSIScore (MSIsensor2, 24Q4) | 1,311 (76%) | 83 have > 20 (DepMap's MSI threshold) |
| CMP msi_status / ploidy (Sanger, 2021) | 610 / 618 | concordance with DepMap MSIScore > 20: 47/47 CMP-MSI have > 20; 7/563 CMP-MSS have > 20; ploidy r = 0.52 (n = 477) |
| any ploidy / any MSI | 1,178 / 1,311 | |
| derived arm-level CN | 1,724 (100%) | §7 |
| purity | 0 | not published; assume ≈ 1 for cell lines (**A5**) |

The 24Q4 ploidy distribution has median 3.06 (IQR 2.33–3.69), and the WGD rate is 72.5% among covered lines.
**Q3:** ploidy/WGD are missing for 40% of matched models (all organoids and post-24Q4 lines).
The options are: (a) restrict to covered models, (b) add a missingness indicator, or (c) use relative-CN-derived proxies (arm CN, CIN-like fraction of genome altered).
Since relative CN is already ploidy-normalised, ploidy matters mainly for expression (dosage) effects.

## 10. Benchmark: PIE X-Atlas/Orion DE tables

**Definitions (A6–A9):**
- **A6, tested:** a pair counts as tested when its row is present.
- **A7, filter:** a feature is present for a target iff its mean is ≥ 0.05 (CP10k ≈ 5 CPM) in the control *or* the perturbed group.
  Features with control mean ≥ 0.05 are tested for **every** target ("consistently tested": HCT116 11,100, HEK293T 12,247). The rest (control mean < 0.05; 9,248 in HCT116) appear only when a perturbed group passes by chance.
  **Q4:** "target measured" in the universe should use the consistently tested features (HCT116 10,368 and HEK293T 11,171 of them are DepMap-expressed); they hold 88% (HCT116) and 94% (HEK293T) of the trans hits.
- **A8, BH scope:** BH is applied *within each target* (recomputed BH matches p_adj to 2e-16).
- **A9, trans:** both genes are placed and on different chromosomes. It covers 94.4% of tested pairs; cis is 5.3%.

| | HCT116 | HEK293T |
|---|---:|---:|
| DE rows (tested pairs) | 102,228,672 | 53,646,719 |
| targets / features | 8,987 / 20,348 | 4,297 / 18,518 |
| features per target (median, range) | 11,355 (11,190–12,679) | 12,468 (12,325–13,134) |
| perturbed cells per target (median, 5–95%) | 144 (45–378) | 190 (52–510) |
| control cells | 165,777 | 218,838 |
| targets whose own gene is tested (expressed) | **3,838 (43%)** | **436 (10%)** |
| on-target LFC < 0 / p_adj < 0.05 & LFC < 0 / p_adj < 0.01 & LFC < 0 | 100% / 89.8% / 86.3% | 100% / 82.6% / 79.1% |
| on-target LFC median (5–95%) | −27.5 (−29.6 to −3.0) | −26.8 (−29.1 to −2.2) |

- **LFC artefact:** in 83.5% (HCT116) and 80.5% (HEK293T) of expressed targets the perturbed-cell mean is exactly 0, and LFC then reflects only a tiny pseudocount (−25 to −30).
  |LFC| > 10 occurs in only 0.017% of all rows, but in 1.3% of HCT116 p_adj < 0.05 hits.
  **Use sign and p (or Ueffect), not LFC magnitude, and cap/winsorise LFC.**
- **Targets of non-expressed genes** (own gene not tested; their DepMap expression is also ≈ 0, median log2(TPM+1) ≈ 0.02) are natural negative controls.

**Significant non-self hits (per target; "trans" excludes cis and unplaced):**

| | HCT116 | HEK293T |
|---|---:|---:|
| trans pairs p_adj < 0.05 (targets with ≥ 1 / ≥ 10 / ≥ 100) | 224,024 (5,831 / 1,196 / 299) | 88,233 (2,244 / 531 / 172) |
| trans pairs p_adj < 0.01 (≥ 1 / ≥ 10 / ≥ 100) | 121,094 (4,811 / 730 / 206) | 30,425 (1,628 / 235 / 71) |
| per-target trans hits at 0.05: median / 90th / 99th pct / max | 1 / 13 / 687 / 5,285 | 1 / 14 / 490 / 4,319 |
| share of trans hits from top 1% / 10% of targets | 55% / 92% | 48% / 94% |
| cis pairs p_adj < 0.05 | 12,900 | 5,484 |
| trans p_adj < 0.05 with feature in DepMap expression | 198,732 | 80,854 |
| fraction up-regulated (trans, 0.05) | 47.8% | 57.0% |
| mean trans hits (0.05): own gene expressed vs not | **51.0 vs 5.5** | **34.4 vs 19.0** |
| share of all trans hits from non-expressed targets | 12.7% | **83.0%** |
| Spearman(target_ncells, n trans hits) | −0.30 | −0.18 |
| Spearman(on-target LFC, n trans hits) | −0.02 | −0.03 |

- **Overlap of targets:** 3,277 targets are shared, 5,710 are HCT116-only and 1,020 HEK293T-only (10,007 in the union).
  Across shared targets the Spearman of on-target LFC is 0.48, and of the number of trans hits 0.35.
- **Universe preview** (not the final definition): targets linked to DepMap CN × consistently tested trans features in DepMap expression, excluding self.

  | target filter | HCT116 targets / pairs / positives (0.05) / rate / % down | HEK293T |
  |---|---|---|
  | in DepMap CN | 8,727 / 85.7 M / 182,184 / 0.21% / 61% | 4,074 / 43.1 M / 75,328 / 0.17% / 47% |
  | ≥ 20 lines with \|log2 CN\| > 1 | 6,357 / 62.4 M / 131,245 / 0.21% / 61% | 3,041 / 32.2 M / 58,409 / 0.18% / 46% |
  | in DepMap CN & own gene expressed | 3,802 / 37.3 M / 162,397 / 0.44% / 61% | 431 / 4.6 M / 12,481 / 0.27% / 61% |
  | ≥ 20 lines \|log2\|>1 & own gene expressed | 2,666 / 26.2 M / 116,560 / 0.45% / 60% | 314 / 3.3 M / 9,548 / 0.29% / 62% |

## 11. Things that are surprising or concerning

1. **HEK293T benchmark is dominated by targets that are not expressed.** 90% of HEK293T targets are not expressed in HEK293T, yet they produce 83% of its trans hits (mean 19 per target vs 34 for expressed targets).
   In HCT116, non-expressed targets average 5.5 hits (vs 51), a plausible null level. In HEK293T, hits from non-expressed targets have a different up/down balance (59% up).
   This suggests off-target, guide or batch effects, or cell-number effects. **Q5:** restrict regulators to targets whose own gene is expressed (and, optionally, has significant knockdown), and use non-expressed targets as an empirical null.
2. **HEK293T is not a DepMap model** (HEK293 is). HCT116 is.
3. **Trans hits are concentrated:** 10% of targets give > 90% of hits, and a few (max > 5,000 hits) are global-stress perturbations. Pair-level enrichment will be driven by these few targets unless it is stratified or weighted per target.
4. **LFC magnitudes are not usable for zero-mean groups** (pseudocount artefact).
5. **Ploidy/WGD/MSI are frozen at 24Q4** and missing for all organoids. Purity is not available.
6. **"Variable CN" needs a strict definition** (§6), because aneuploidy makes nearly all genes vary.
7. **Release not md5-verified** against the portal (bot protection); the evidence is strong but indirect.
