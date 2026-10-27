# REPORT.md independent verification (fresh agent, 2026-10-27)

The verifier could not write this file itself (the harness blocks subagents from writing report files), so the
PI saved this condensed record of its findings. All listed edits were applied to REPORT.md in the next commit.

**Verified:** roughly 120 numbers in REPORT.md match `results/final/` (summary.md, the per-run `metrics` JSONs),
`results/tuning/`, `results/m4/` and `results/tuning_log.tsv`.

**Errors (fixed):**
1. The test count was 185; 187 tests are collected.
2. PureCN missingness was given as 40%; in S2 it is 30% (437/1,436).
3. "M1 1.59" was unlabelled and is a TRAIN value.
4. The "logged override" for the LOCK does not exist; final_eval.py refuses to run while the LOCK exists.

**Over-claims (fixed):**
- **Sign concordance.** The p-value comes from a within-regulator null with mean 0.647, not 0.541, so the excess
  is about 7 points. The result is not significant in S1 (p = 0.13) or in the no-split-sibling subset
  (p = 0.48).
- **"Replicates".** The pre-registered permutation p is 0.068, and the paired gain over NC1 is not significant
  (CI −0.66 to 2.85).
- **"No within-regulator information".** Rescoped to "no *detectable* target-ranking information": the C1 CI is
  wide, and the CMH bounds any lift at about 15%.

**Imprecise or under-reported (fixed):**
- Runtime and memory figures, and the evaluation count (55).
- Call-set size (10.7 M genome-wide) and the fact that it is identical across the G1/A3 variants.
- Sign of the leave-HCT116-out Δ.
- B4/B4e/B2 exceed the CN-only OR on TEST.
- The HEK293T context layer is not significant.
- The no-split-sibling subset is weaker.
- No multiplicity adjustment across secondary TEST metrics.
- Reproduction commands: PIE input, paralog/class builds, runtime table, and the warning not to rebuild after
  the freeze.
- Per-method runtime and sign-concordant OR columns.
