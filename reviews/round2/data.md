# Round 2 review: data (DESIGN v2 §2–§4, §6, §11)

I re-derived every number below. To do so I rebuilt S2 from §3: matched cell lines, minus `EngineeredModel` lines, keeping one model per `PatientID` and preferring the one with PureCN.

## Round-1 items

- **Handled correctly:** #1–#6, E1 (task T2), E3, E4, E5.
- **Not dispositioned:**
  - #8, the X arms in LOCO. Xp–Xq r = 0.63, which is sex-driven.
  - #9 and #10. These are optional.

## Checks

- **S2.**
  - `EngineeredModel` exists. Its values: drug adapted 8, engineered clone 3, knockout 1, knock-down 1, knock-in 1. All other lines are NaN.
  - Every drug-adapted line and A549_CRAF_KD has a parental line in the panel. A name screen finds no unflagged derivatives.
  - S2 = **1,436**: 14 engineered lines and 59 same-patient duplicates are dropped. HCT116 is in S2.
  - Coverage: ploidy 999; MSIScore 1,255 (CMP adds 0); `Stranded` 1,436 (100%).
- **Frozen flags.**
  - There are 11,425 autosomal variable regulators, insensitive to the clip.
  - On X, male lines sit at log2 ≈ −1. 320 lines per X gene pass the threshold on raw CN vs 58 after sex-centring. The flag count (716/756) agrees only by coincidence of the threshold.
  - Requiring CN costs 158 of 10,368 HCT116 targets. This is fine.
- **Arm-triplet split.** It is feasible.
  - 42 arms hold HCT116 universe regulators (5,623 in total, 1–431 per arm).
  - Across 2,000 simulated draws, TEST holds 32–35% of regulators but 26–40% of positives (90% interval). The top 1% of regulators carry 58% of positives; PHB1 alone has 4,099.
- **TP53.** The 26Q3 mutation files are portal-only. 24Q4 figshare has `OmicsSomaticMutations`, `…MatrixDamaging` and `…MatrixHotspot`, all md5-checkable. They cover 1,255 of 1,436 lines; the missing lines are mostly CNS/Brain (80 of 181).
- **Platform.** The RNA `SequencingPlatform` field is NaN for 93% of lines. Its only level with ≥ 15 lines has 67 lines.
- **Lineage rule.** Only some heterogeneous lineages have subtypes with ≥ 15 lines:

  | lineage | split | residual |
  |---|---|---|
  | Lymphoid | Mature B 65, B-ALL 21, T-LL 15 | 28 |
  | Myeloid | AML 40 | 13 |
  | CNS/Brain | glioma 118, embryonal 21 | 10 |
  | Bone | osteosarcoma 44, Ewing 24 | 12 |
  | Soft Tissue | none | – |

  "Parent or other" is ambiguous for residuals under 15.

## Changes

1. **MUST: name the TP53 source and status rule.** Use **24Q4** (DOI 10.25452/figshare.plus.27993248) and record it in DATA.md. TP53-mutant = Damaging > 0 OR Hotspot > 0. Note that the missingness indicator mostly marks the new CNS lines.
2. **MUST: specify the X rule.** Flag variable X CN on sex-centred log2 CN, excluding unknown-sex lines. Also sex-centre Xp/Xq in LOCO, or drop them.
3. **SHOULD: drop `SequencingPlatform` as a covariate.** Keep `Stranded`.
4. **SHOULD: make the lineage residual rule explicit.** Residuals become a "<lineage>-other" level whatever their size; with 10–28 lines each is estimable. Lineages with fewer than 15 lines go to one "other" level. It holds 47 lines, NaN lineage included. Add explicit levels for GrowthPattern "Unknown" (8) and media NA (38). Freeze the level table on P.
5. **SHOULD: log the exact n.** Replace "≈ 1,450" with the exact frozen n. Run T3 before the split, because the arm ordering depends on the frozen flags.
6. **COULD: improve ploidy imputation.** Lineage-median imputation overstates ploidy for the quiet and MSI-H lines that lack PureCN. Instead, impute from MSI plus genome-wide arm burden, fitted on the 999 lines with ploidy.
7. **COULD: two small additions.** A sensitivity excluding the ~20 non-cancerous lines, and keeping Xp and Xq in the same split.

**No change needed:**
- the engineered exclusion and patient dedup;
- the homozygous-deletion indicator at relCN < 0.25 (81% of values below 0.25 are below 0.05);
- the arm-triplet design;
- excluding Y and PAR genes.
