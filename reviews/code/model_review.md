# Independent review: model engine (benchmark-blind)

Reviewer: independent agent, 2026-10-27. Scope: DESIGN v3.1 §1, §3–§5, §7, §8.6; RESEARCH_TREE D3;
docs/INTERFACES.md; docs/MODEL.md; `src/depmaptr/{io,covariates,assoc,gates,attribution,simulate,pipeline,baselines,calibrate}.py`,
`scripts/run_model.py`, `results/variants.yaml`, `tests/test_model_*.py`, `tests/test_sim.py`, and the real score files in
`data/results/scores/`. I did not open any benchmark file (bench/, pie/, pie_x_atlas_orion).

All checks use my own code in the session scratchpad
(`/private/tmp/claude-503/-Users-colmclaude-ClaudeExperiments-DepMapTR/130bbdd0-af4a-481e-a7c7-41816a9e851c/scratchpad/`:
`indep.py` rebuilds the panel, Z, CN transform, arms, burden and D_c from the raw processed files without importing `depmaptr`;
`chk1.py`–`chk9.py` hold the comparisons; `ablate.py` holds the SIM ablations). I changed no repo file apart from this one.

## Verdict

The core engine (M1 t, G1 cis t, s_A, G2, the gated composition, the baselines and NC1) reproduces exactly from independent
code. One **blocker** remains: the A3 run with co-expression weight w_co > 0 has NaN scores for 28% of regulators. Fix it
before any tier-4 run is evaluated. λ_GC ≈ 2.1 shows no sign of a bug.

---

## BLOCKER

### B-1. A3 with w_co > 0: NaN α for every multi-gene segment that contains a gene without expression
- **Evidence (real run `a3_r095_s2_w05_g2_thr1_w0`):**
  - 63,152,238 NaN scores (28.5% of all pairs); `t` has only the usual 11,716.
  - 3,455 of 12,141 regulators have an **all-NaN** α row; no row is partially NaN.
  - 2,039 of those regulators are expressed in HCT116. The harness will rank every one of their pairs at the minimum.
  - Per chromosome: chr1 347/1,568 rows, chr8 198/640, chr17 237/983, chrX 231/716.
- **Mechanism (`pipeline.run`, `attribution.segment_alpha`):**
  - `dfco_c` is initialised to NaN and filled only for genes with an expression column. CN-only genes (426 regulators, plus
    non-regulator CN genes) keep `df_co = NaN`.
  - `log_abf(z=0, df=NaN)` then returns NaN, because λ = σ₀²/(σ₀² + 1/NaN). That NaN enters `ll`, and the segment's
    `logsumexp` makes α NaN for every member and every target.
- **Independent confirmation (`chk5.py`):**
  - I rebuilt chr8, chr18 and chr21 from scratch: LOCO residualisation, greedy segmentation, t_gB for all CN genes, cis t,
    z_co, π and α. My segment counts equal the run's `n_segments_by_chrom` (122, 57, 87).
  - On finite entries, α matches the stored α to 3e-8.
  - The NaN pattern is reproduced exactly by mimicking `df_co = NaN`. Every all-NaN row lies in a multi-gene segment that
    contains a CN-only gene.
- **Why the tests missed it:**
  - SIM gives every gene an expression column.
  - `test_cosegment_attribution` uses w_co = 0.
  - `test_real_score_files` covers only M1, B0–B6 and NC1, not variant runs. Its finite-fraction threshold (> 0.95) would have
    caught this.
  - MODEL.md §9 reports the run with no NaN note.
- **Impact:** tier-4 options `a3_r095_s2_w025` and `a3_r095_s2_w05` would be scored as if a quarter of the regulators had
  been dropped. This would bias the adoption decision against them, and it would spend part of the 30-configuration cap
  on a broken score.
- **Fix:** when a gene has no expression, the co-expression factor should be absent (ABF^w = 1). Set its `ll` contribution
  to 0, e.g. `np.where(np.isfinite(df_co), w_co*log_abf(...), 0)`. Then:
  - add a SIM with some CN-only genes;
  - extend `test_real_score_files` to every run directory, with finite share ≥ (1 − (A=B + CN-only-B2)/R·T).

  On the finite entries my spec-faithful α equals the stored α, so the fix should change nothing else.

---

## MAJOR

### M-1. A3 + G2: G2 is not enforced for 1,384 non-expressed regulators (DESIGN-literal, undocumented consequence)
- DESIGN §5.6 lets w_G2 enter only through π_g. Two §5.4 rules then switch G2 off:
  - singleton segments get α = 1;
  - when Σπ = 0, π falls back to uniform.
- In `a3_r095_s2_w05_g2_thr1_w0` (G2 w = 0), 5,061 regulators are not expressed in HCT116:
  - 1,416 are all-NaN (B-1);
  - 2,261 are correctly zeroed;
  - **1,384 keep non-zero scores**. Of these, 1,020 are singletons with α ≡ 1, so they get their full M1 t. The rest sit in
    segments where every member is non-expressed, so π is uniform.
- Example: on chr21, 32 of the 41 non-expressed regulators are singletons.
- I confirmed that score == t·α exactly, with no G2 multiplier, as the code intends.
- **Consequence:** "A3 + G2" vs "G2 only" mixes the attribution effect with *removing* the G2 penalty from about 27% of
  non-expressed regulators. That runs against the stated rationale for G2 (a non-expressed regulator cannot be knocked down;
  §5.3).
- **Recommendation:** a PI decision note, either:
  - apply w_G2 as a multiplier for singletons and for the Σπ = 0 fallback; or
  - keep the literal rule and document it in MODEL.md §6.

### M-2. SIM does not exercise two of its listed confounders; one adjustment test is fragile
I removed each adjustment in turn by monkeypatching `ablate.py` (repo untouched), then reran the SIM tests on seeds 1–3.
Each cell shows the result per seed (P = pass, F = fail):

| adjustment removed | arm co-selection | lineage | sign rule | AUROC | A3 | null λ |
|---|---|---|---|---|---|---|
| none | PPP | PPP | PPP | PPP | PPP | PPP |
| lineage in Z | PPP | **FFF** | PPP | PPP | PPP | PPP |
| LOCO arm block | **FPF** | PPP | PPP | PPP | PPP | PPP |
| W_B entirely | **FPF** | PPP | PPP | PPP | PPP | PPP |
| deletion indicator only | PPP | PPP | PPP | PPP | PPP | PPP |
| sex in Z | PPP | PPP | PPP | PPP | PPP | PPP |
| s_A rule (s = +1) | PPP | PPP | **FFF** | PPP | PPP | PPP |

- **Homozygous deletions** (listed in §7) are simulated but inert:
  - no planted regulator's CN correlates with a target's deletion;
  - no test asserts anything about deletion targets.

  So removing the §4 deletion indicator goes undetected.
- **Arm co-selection:** with the LOCO arm block removed, the test still passes on seed 2. The burden covariate partly absorbs
  the co-selection, and the assertion has little margin. Two remedies:
  - use more seeds, or assert a contrast (fp_M1 − fp_noarms);
  - add a monkeypatched negative control to the suite.
- **Not in SIM:**
  - planted X regulators or any sex-dependent confounding (removing sex passes);
  - genes without expression (which hid B-1);
  - burden, TP53 or MSI confounding.

  §7 does not require these, but nothing tests the burden, TP53, MSI or sex adjustments at all.
- The A3 test uses only w_co = 0, and it has no CN-only or non-expressed members.
- What works: the lineage and sign-rule tests are real tests (they fail when the adjustment is removed), and null calibration
  holds.

---

## MINOR

### m-1. Hinge variant on X regulators
`group_t` splits `x_all` at 0 without sex-centring. Male X (log2 ≈ −1) therefore falls entirely in the "loss" part, and
male gains up to log2 0 count as loss-slope variation. The sex dummies absorb the mean, but not the hinge location. DESIGN
says "hinge at copy-neutral". For male X, copy-neutral is −1. This affects 716 X regulators under `hinge` only, and it is
undocumented. Fix: either hinge X regulators on sex-centred x, or document the choice.

### m-2. Lineage level count in MODEL.md
- There are **32** lineage levels: 31 dummies plus the reference Lung. MODEL.md §3 and §7 say "31 levels" and "the 31 Z
  lineage levels".
- The NC1 strata are the 32 levels.
- My independent rule gives identical labels.

### m-3. B0 meta seed field
When the default seed is used, `meta.json` top-level `seed` is `null`. The seed actually used (20261027) appears only in
`diagnostics.seed`.

### m-4. B3 sign
B3 is ranked by |cis t × SD(x_A)|, so negative-cis regulators rank as hubs. This is negligible on S2: the most negative
cis t is −3.64, against ≈ 20 for typical regulators.

### m-5. NC1 same-gene entries
NC1 has 0 NaN scores: the 11,713 A = B entries are finite after permutation. This is harmless because the harness masks
same-chromosome pairs, but it differs from M1's NaN convention.

### m-6. §5.2 empirical null and call set not implemented
Efron z′, BH(z′) ≤ 0.05 and BY are not implemented. MODEL.md flags this as deferrable (§12). Note that with λ ≈ 2 any
nominal-t call set would be meaningless until it exists.

### m-7. Documented deviations (agree they are documented)
- RINT arm CN keeps log2.
- Remainders of split parent lineages below 15 lines (Myeloid 13, Bone 12, CNS/Brain 10) stay as their own level.
- Target standardisation excludes same-chromosome entries.
- A3 segments use LOCO D_c under all arm handlings.
- In loao, a regulator on an arm with < 20 genes keeps all arms, including the sibling arm.

---

## Verified (no issue found)

### 1. Covariates (S2), rebuilt independently (`chk1.py`)
- Lineage rule: identical labels for all 1,436 lines (32 levels).
- Z: my 58 columns (before dropping duplicates) span the same space as the repo's 56. The residual of each on the other is
  ≤ 2.7e-13, and the rank of [1, Z] is 57 both ways.
- Encodings:
  - missing indicators;
  - ploidy/WGD: lineage-median imputation;
  - MSI > 20; TP53 altered;
  - sex: Female and Unknown-or-NaN;
  - media base levels ≥ 15, with the rest and missing pooled as `other`, reference RPMI;
  - serum-free flag and its missing indicator;
  - growth pattern, reference Adherent;
  - Stranded.
- Platform is dropped (D3). `wgd_missing` ≡ `ploidy_missing` and `tp53_missing` ≡ `msi_missing` are dropped as exact
  duplicates.
- In-scope CN genes: the same 18,998.
- Arm CN: the same 41 arms. Arm values are the median of the active log2-clip transform over arm genes, and Xp/Xq are
  sex-centred within Male / Female / Unknown-or-NaN. Max abs difference is 0.
- D_c for chr1, chr8 and X spans the same space as my [1, Z, LOCO burden (autosomal arms ∉ c), arms ∉ c], with residual
  ≤ 2.6e-13.

### 2. M1 t vs my statsmodels OLS (`chk2.py`, `chk3.py`)
I fitted y_B ~ D_c + x_B + 1(relCN_B < 0.25 if ≥ 3 lines) + x_A from raw data on 17 real pairs. They include:
- X regulators: STS→RCC2, EDA→CDKN2A;
- homozygous-deletion targets: CDKN2A (426 lines), PTEN, TP53, and 5 random deletion targets;
- the negative-cis regulators SLAMF1, LRP1B, DEFB116;
- a CN-only regulator;
- MYC, ERBB2, CCNE1, MDM2, KRAS.

Results:
- Max |Δt| = 5.5e-8, which is float32 precision.
- df = 1,335 with the deletion indicator and 1,336 without, matching n − rank(D) − k_B − 1.

### 3. G1 and the s_A rule
- G1 cis t (y_A ~ x_A + D_c) matches stored `cis_t` on 16 regulators (≤ 1e-4).
- Exactly the 3 regulators with t_cis ≤ −3 get s = −1. The most negative t among s = +1 regulators is −2.46.
- CN-only regulators have t_cis = 0 and s = +1.

### 4. G2
- Model.csv lists ACH-000971 as "HCT 116" and ACH-001085 as "HEK-293".
- Context expression is read from the full matrix; 7,080 regulators have log2(TPM+1) ≥ 1.

### 5. `g1_t30_w0_g2_thr1_w0` composition
- `t` and `cis_t` are bit-identical to M1.
- score == t·1(t_cis > 0 ∧ t_cis² ≥ 30)·1(expr ≥ 1), with 0 max difference. 7,457 regulators pass G1 and 6,481 pass both.

### 6. A3 formula (excluding B-1 entries)
- λ = σ₀²/(σ₀² + 1/df) with σ₀ = 0.05; ABF^w_co with w = 0.5.
- z_co = atanh(r)·√(df − 1), with r = partial corr(y_g, y_B | x_g, W_B, D_c).
- π ∝ exp(t_cis/2)·w_G2, with the uniform fallback.
- Singleton segments get α = 1.
- α is normalised over **all** CN genes in the segment (chr21: 205 genes, 56 of them regulators).
- My α matches the stored α to 3e-8. G2 is not double-applied (score = t·α).

### 7. Baselines and NC1 (`chk7.py`, `chk8.py`)

| control | definition checked | match |
|---|---|---|
| B0 | U(0,1) | yes |
| B1 | Pearson t, df n − 2; naive cis for s_A | ≤ 1e-6 |
| B2 | y_A vs y_B; CN-only rows NaN; s = +1 | ≤ 1e-6 |
| B3 | M1 cis × SD(x_A), constant per row | ≤ 2e-6 |
| B4 | SD(y_B) | ≤ 2e-7 |
| B4e | HCT116 y_B | exact |
| B5 | rank(\|B3\|)·rank(B4e) | rel ≤ 6e-8 |
| B6 | own-arm CN t, LOCO D_c + W_B, checked with statsmodels | 5 d.p. |
| NC1 | one SeedSequence(20261027).spawn(1)[0] permutation within the 32 lineage strata, applied to the D_c-residualised x_A; Ỹ and W̃_B fixed | 11 pairs, 5 d.p. |

### 8. λ_GC ≈ 2.1 (`chk9.py`): structure, not a bug
- The engine is exact against OLS, and NC1 gives λ = 0.96, so df and residualisation are correct.
- Inflation is broad: the 10th–90th percentile of per-regulator mean t² is 1.5–3.2. No single bad block drives it.
- Each M1 regulator row correlates at median r = 0.85 with the B6 row for its own arm.
- On the same 400 rows, λ is 1.23 for focal and 2.49 for B6.

So most of the inflation is A's own-arm dosage acting on genome-wide expression: the §1 LOCO caveat. The residual 1.2
under focal is plausibly genuine polygenic trans signal plus unmodelled confounding. It is not an engine artefact. The
§5.2 empirical null is the intended remedy (m-6).
