# Round 1 review: data foundation and data-level design choices

Reviewer scope: independent audit of `docs/DATA_PROFILE.md` claims against `data/processed/`, plus the data choices in DESIGN v1 §2, §3, §4, §6. All checks were re-derived with throwaway scripts outside the repo, using `.venv/bin/python` on the processed matrices and PIE caches. No repo files were edited apart from this one.

## (a) Audit findings

### Errors and omissions

**E1. The PureCN provenance numbers are wrong (DATA_PROFILE §9).** The profile says "905 from the *same* WGS profile, the rest model-level". Of the 1,037 lines that have PureCN values:
- 687 are `same_profile` (WGS);
- 350 are `model_level_24Q4`, and 349 of those are WES.

905 is the number of matched cell lines whose 26Q3 CN ProfileID exists in 24Q4 `OmicsSignaturesProfile`. 218 of those 905 have MSIScore but NaN Ploidy/WGD/CIN/LoH/Aneuploidy at both the profile and the model level. So there is no fallback bug in `build_depmap_tables.py` (I checked), but the sentence is wrong. The 687 in §4 is correct.

**E2. Omission with design consequences: HCT116 is not in S1.** ACH-000971 has MSIScore 89.1 but no PureCN ploidy/WGD/CIN/LoH in 24Q4, at either the profile or the model level. Missing PureCN is not random:
- Among lines with an MSIScore, PureCN is missing for 37/83 MSI-H lines (45%) but only 237/1,228 MSS lines (19%).
- Lines outside S1 have quieter genomes: the median fraction of arms with |log2| > 0.3 is 0.27 outside S1 vs 0.37 in S1.
- 97 of the 154 CNS/Brain lines are missing.

So S1 under-represents near-diploid MSI-H, HCT116-like lines, and it drops the benchmark line itself. DESIGN §4 cites "HCT116 itself is MSI-H" to justify the MSI covariate, and §5.4 proposes weighting toward HCT116-like lines. Both sit badly with a sample set that excludes these lines.

**E3. The variable-CN counts are not those of the sample set the design uses.** DATA_PROFILE §6 and §10 compute variability on all 1,724 models, organoids included. DESIGN §3 computes it "in the sample set". My recomputation (autosomal genes, ≥ 20 lines with |log2 relCN| ≥ 1):

| sample | variable genes |
|---|---:|
| 1,724 matched models | 12,999 (matches the profile) |
| 1,509 cell lines | 11,748 |
| S1 (1,037) | 7,796 |

So for S1, the universe preview overstates the regulator set by about 40%. More importantly, a fixed count threshold depends on n. If S1 vs S2 stays a tunable, **the benchmark universe changes with the tuning choice**, and ORs computed on different universes cannot be compared.

**E4 (minor). PAR genes link to the Y copy.** IL9R and WASH6P link to `"ENSG…_PAR_Y;SYMBOL (id)"`, and VAMP7 links only to the `_PAR_Y` column. `profile_benchmark.var_flag` takes `split(";")[0]`, which is the Y copy. `test_pie_map_one_to_one` passes because it tests the joined strings. X regulators are excluded, so this has no effect on regulators, but it would give these targets the wrong cis CN (x_B).

**E5 (minor). The CN floor is inconsistent.** DESIGN §3 uses 0.05, while the code and the arm CN use 2^-5. The two affect about the same share of values (0.086% vs 0.085% fall below them), but one constant should be used.

### Claims confirmed

- **Matching:** 1,724 matched models (1,509 cell lines, 215 organoids). S1 = 1,037 cell lines, all with MSIScore. WGD = 752/1,037 (72.5%), consistent with ploidy (mean 2.14 without WGD, 3.46 with). 14 S1 lines have CN and RNA from different ModelConditions.
- **CN scale:** linear relative CN. Min 9.5e-7, max 260. ERBB2, MYC and CDKN2A controls reproduce.
- **24Q4 vs 26Q3 CN agreement (new check):** the per-line log–log Pearson of 26Q3 relative CN vs 24Q4 PureCN absolute CN has median 0.92 for both same-profile WGS lines (n = 686) and model-level WES lines (n = 350). The 37 lines with r < 0.5 are flat diploid genomes (SD of log2 CN 0.10 vs 0.40 overall), not sample mismatches. **The 24Q4/26Q3 mismatch risk is low.**
- **Gene mapping:**
  - 99.3% of HCT116 targets linked to CN have the identical symbol in DepMap. The other 58 are HGNC renames and resolve correctly through Entrez (e.g. SLC22A18 → "SLC67A1 (5002)", NDUFA4 → "COXFA4 (4697)").
  - The chromosome assignments of the PIE map and the DepMap annotation agree for all 18,222 expression links. The 3 CN-link mismatches are E4.
  - Cross-dataset identity check: PIE HCT116 control mean vs DepMap HCT116 log2(TPM+1), Spearman **0.89**. Only 0.6% of consistently tested features have DepMap HCT116 expression below 0.5.
- **"Consistently tested":** exactly equivalent to `ref_mean ≥ 0.05`, with 0 discordant features (HCT116 11,100; HEK293T 12,247).
- **"Own gene expressed":** own-feature-tested is exactly equivalent to the own feature's `ref_mean ≥ 0.05`, with 0 discordant (3,838 and 436 targets). This is a clean definition.
- **Trans/cis:** landmark genes are placed correctly. PAR genes count as X (A2).
- **Arm CN:** the recomputed HCT116 8q (0.5908) matches the parquet exactly. HCT116 shows +8q and +17q, and male X at −0.97, consistent with its karyotype. The pan-cancer inter-chromosomal arm |r| is ≤ 0.26 in S1, but within-lineage correlations are stronger (Bowel 18q–20q r = −0.46, n = 102). Lineage effects plus LOCO arms are therefore needed, as designed.

## (b) Design changes

**1. MUST: freeze the universe flags on a fixed reference panel.** Compute the regulator variable-CN flag and the DepMap "measured" flag once, on a panel that no tunable can change. I recommend the deduplicated cell-line panel (see #3). Commit them next to `regulator_split.tsv` before any tuning. Sample-set choices (S1/S2) then change only the model, never the universe. Rationale: E3.

**2. MUST: make S2 (all cell lines) the primary sample set and S1 a sensitivity analysis.** Rationale: E2. S1 drops HCT116, 45% of MSI-H lines and 31% of power (t-statistics scale as √n, a factor of about 0.83). The brief's "~1,000 with both" refers to lines with matched CN and RNA, not to PureCN coverage, so it does not justify S1. Handle covariates in S2 as follows:
- **Ploidy, WGD, CIN, aneuploidy, LoH:** add a missingness indicator plus an imputation. Fit it in S1 from relative-CN features, which carry real information: the fraction of altered arms has Spearman 0.62–0.64 with CIN, aneuploidy and WGD. Candidate features are the arm-FGA, the spacing of segment levels and the per-line SD.
- **MSI:** use MSIScore, else Cell Model Passports `msi_status`, else an indicator.

If the PI keeps S1 vs S2 as a tunable, #1 still applies and both results must be reported.

**3. SHOULD: remove pseudo-replicates.** Keep one model per PatientID: 58 patients contribute 70 extra cell-line models (35 extra in S1). Also drop drug-adapted and engineered derivatives when the parental line is present, e.g. A-375 ×3, HT-144 ×3, RPE1-ss ×3, HCC-827-GR5, A549_CRAF_KD. Rationale: near-identical genomes inflate the evidence for their shared CN events. They also break exchangeability in the NC1 within-lineage permutation and in the bootstrap.

**4. SHOULD: coarsen the held-out split blocks.** CN correlation within an arm decays slowly. The median gene–gene r(CN) in S1, across chromosomes 1, 3, 5, 8, 12 and 17:

| distance | median r(CN) |
|---|---|
| 5–10 Mb | 0.56–0.82 |
| 10–20 Mb | 0.44–0.71 |

With 5 Mb bins, most TEST regulators have a TRAIN neighbour with a nearly identical CN profile, so DepMap-side scores are shared across the split. Proposed changes:
- Use contiguous blocks of ≥ 25 Mb that never cross an arm boundary.
- Report the distribution of max |r(CN)| between each TEST regulator and any TRAIN regulator.
- Hub concentration matters: the top 1% of regulators give 55% of HCT116 trans hits. Pre-declare one balancing rule and run it once, logged. For example, sort blocks within each chromosome by PIE positive count and assign them systematically. Also pre-declare a hub-excluded sensitivity metric (regulators above the 99th percentile of hits removed). Re-drawing seeds after looking at labels is not allowed.

**5. SHOULD: add RNA strandedness as a technical covariate.** 44% of S1 lines are stranded vs 66% of non-S1 lines, so library protocol is confounded with sample set and onboarding era. It shifts quantification of genes with antisense overlap. If it is ignored, protocol batch × lineage composition can create shared target programmes.

**6. SHOULD: fix CN leverage at the floor.**
- Use one floor constant (E5).
- With a floor of 0.05, homozygous deletions sit at log2 −4.3, about 12 median per-gene SDs below 0. 0.5/99.5% winsorisation removes them only for genes with homozygous-deletion frequency < 0.5%. Passenger homozygous deletions at fragile sites (frequency 0.5–5%) will then dominate OLS slopes.
- Pre-declare a floor of 0.25 (log2 −2) as a variant next to the rank-INT variant. Report the maximum leverage per regulator.

**7. SHOULD: make the covariate encodings explicit and complete.**
- **GrowthPattern:** also has Spheroid (76 cell lines, mostly CNS stem-media lines outside S1) and Unknown (9).
- **Sex:** "Unknown" (70 cell lines) becomes its own level.
- **Lineage:** NaN (8 lines, including HEK293) goes to "other".
- **Media:** parse the base medium from `mc_FormulationID`. Approximate levels among cell lines: RPMI 725, DMEM 216, MEM 78, IMDM 54, EMEM 51, F12 50, McCoy's 26, AlphaMEM 23, L-15 14. Add DMEM:F12, a stem/serum-free flag, "other" and NA (41). Optionally add serum % bins. Compute the merge thresholds on the fixed panel.

**8. COULD: X arms and PAR genes.**
- Xp/Xq arm CN in the LOCO set is almost a sex dummy (male lines ≈ −1). Either sex-centre the X arms or drop them from LOCO.
- Fix E4 by preferring the X CN column for PAR genes.

**9. COULD: run a sensitivity analysis excluding the 36 models (14 in S1) whose CN and RNA come from different ModelConditions.**

**10. COULD: stratify sign concordance by how a regulator's CN varies.** In S1, 2,596 variable regulators are gain-only (< 5 lines with log2 ≤ −1) and 792 are loss-only. "CN loss ↔ knockdown" is an extrapolation for gain-only regulators, so report concordance separately for loss-supported regulators.

### No change needed

- **24Q4 PureCN vs 26Q3 CN:** concordant (see audit). Keep the `sig24_source` flag as a covariate-provenance indicator.
- **Variable-CN rule (≥ 20 lines, |log2| ≥ 1):** sensible once #1 fixes the panel.
- **DepMap "measured" = SD > 0:** a no-op (0 of 10,368 HCT116 universe features fail it), but harmless, because the PIE consistently-tested filter already restricts targets to expressed genes. 45 of these features are expressed (log2(TPM+1) ≥ 1) in fewer than 10% of lines; the robust variant covers them.
- **Positives (`p_adj < 0.05`, BH within target), consistently-tested target filter and expressed-regulator secondary universe:** the definitions are verified exactly.
- **HCT116 primary, HEK293T secondary:** agreed. Report HEK293T on the expressed-regulator universe, which has only about 314 regulators.
- **Leakage from HCT116 being a DepMap line:** none. Its CN/RNA profile is DepMap data, not benchmark labels, and it is 1 of about 1,500 lines. Keep it in S2. Optionally run one sensitivity without it, and cap the per-line weight in the HCT116-context-weighted variant so that HCT116 cannot dominate.
- **Excluding organoids:** agreed.
