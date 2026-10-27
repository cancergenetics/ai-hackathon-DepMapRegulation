# Pre-flight review of the one-shot TEST evaluation (independent reviewer, 2026-10-27)

Scope: PROMPT.md; DESIGN §5.5–§8; RESEARCH_TREE D0–D8; `results/final/selected.DRAFT.yaml`;
`scripts/final_eval.py`; `src/depmaptr/{benchmark,controls,report}.py` (the relevant parts);
docs/BENCHMARK.md §10–11; `results/final_mock/`. I did not run final_eval without `--mock` and did not touch
`~/DepMapTR_sealed`. The only repo file I changed is this review.

Independent re-run: I cloned HEAD (2febf06, clean tree) into the scratchpad, symlinked `data/`
read-only sub-dirs, and ran `final_eval.py --mock --selected selected.DRAFT.yaml` there
(`DEPMAPTR_ROOT=<clone>`, under `/usr/bin/time -l`). This touched no repo file.

## Verdict: **GO, conditional.** No blockers. Fix or explicitly accept M1–M3 before the run

Nothing found would crash the run or corrupt it. The selection matches D7/D8. Every run exists and is
consistent. The guards, LOCK and seal logic work, and the mock reproduces bit-for-bit at HEAD.

M1–M3 are rated major because they cannot be fixed after the run: any later fix would need a second read of
the TEST labels. M1 is a cheap guard. M2 and M3 change frozen files, so each needs a decision note,
`record_frozen_change.py`, a re-run of the mock and `tests/test_frozen.py`. If the PI chooses to accept M2
or M3 instead of fixing them, the decision note must say so before the run.

## Evidence summary

**(1) Selection vs D7/D8.** `selected.DRAFT.yaml` implements D7 as written:
- selected = `g1_t100_w0__a3_r095_s2_w0`;
- ablations: minus G1 (`a3_r095_s2_w0`), plus context layer (`…__i1_fuse_g2`), focal arm handling;
- baselines B0–B6 including B4e;
- controls: NC1, C2_g2, PC2 and NC2 (`nc2_run` = selected, which has no G2);
- C1 on the selected run and on the context layer;
- sensitivity runs: S1, noHCT116, hidden_factors, cn_outlier_exclusion, tp53_chr17_drop;
- subsets: paralog_exclusion, class_i and hub_excluded, with no_split_sibling added automatically;
- regulator-class tables and the runtime table appear in every run's output;
- both lines; n_perm = 1000.

The context layer has role "ablation". This is correct, because `roles()` de-duplicates and listing it in
`reported` would change nothing. PC1 and SIM are label-free and not needed on TEST.

All 19 unique score runs exist under `data/results/scores/` and match each other:
- every `scores.npz` holds `score` [12141 × 18276] float32;
- the regulator and target arrays are identical in content and order to the selected run's, with no
  duplicates;
- `s_cis` ∈ {−1, +1} (B2 is all +1, as expected for co-expression);
- the fused run holds `pair_score` and `target_weight`, so its C1 runs in exact mode.

The fresh mock's `data_checksums` match the committed mock for all 51 JSONs, so the score files have not
changed since the committed mock.

**(2) PROMPT coverage.** The run reports:
- direction-agnostic OR at K = 0.01, 0.1 and 1 %, with Fisher p, regulator-label permutation p,
  degree-stratified permutation p and the arm-jackknife CI;
- the sign-concordant OR, concordance vs c₀, the sign 2×2 OR and the within-regulator sign permutation p;
- universe N per line, both lines, and TEST-arm regulators only (HCT116 1,817 regulators and HEK293T 798,
  each on 14 arms; no TEST `own_class` or `degree_decile` values are in the tree).

**Runtime is included in summary.md**, in the `## runtime` section: `eval_wall_s`, `score_wall_s` and
`threads` per role × run × line. `results/runtime_table.tsv` also exists. However, the runtime sits in its
own table, not next to the ORs, and the threads column does not record BLAS threads (m8).

**(3) Mock sanity.** The controls behave as expected:

| run | HCT116 OR | HEK293T OR (jackknife CI) |
|---|---|---|
| B0 | 1.06 | 0.62 (0.35–1.08) |
| NC1 | 0.93 | 0.51 (0.20–1.33) |
| B6 | 0.90 | 0.93 |

- The permutation p values of these nulls lie between 0.35 and 0.95.
- PC2 = 113.
- The selected run scores 4.11 on HCT116 (TRAIN 5.0) and 6.96 on HEK293T.

**NC2 = 0.18 is small-sample and structural, not a bug.**
- HCT116: a = 0 of K = 9,289, where the expected a under independence is about 2.8. The OR is
  Haldane-corrected ((0.5 × 9.28M) / (9289.5 × 2767.5)), with one-sided Fisher p = 1.
- Only 39 of the 965 NC2 regulators have any non-zero score after the G1 gate (w = 0).
- The HCT116 top-K comes from 6 regulators (SETBP1, GSTM1, SOX2, FOXO6, ABHD16B, …), which have 2 positives
  between them in total.
- HEK293T is driven by essentially one regulator: GSTM1 contributes 7,105 of the 7,114 top-K pairs.

In short, NC2 cannot tell much about a G1-gated method (m1).

**(4) LOCK and seal.**

Order of operations (in the real run):
1. Guards run first:
   - an existing LOCK, in the working tree or at HEAD, makes the run refuse;
   - `selected.yaml` must be tracked and clean;
   - there must be no TEST JSON or summary under `results/final/` and no TEST rows in the tuning log;
   - `results/FROZEN.sha256` must exist, with no drift;
   - there are at most 3 ablations, and every run directory must resolve.
2. LOCK is written and committed **before** any sealed file is chmodded or read.
3. Only the 4 TEST label files are set to `u+r`, and each sha256 is checked against the sealed block of
   `FROZEN.sha256`.
4. A `finally` block removes the read bits again (back to 000), writes the summary, marks the LOCK
   COMPLETED or FAILED, and commits.

Checks behind this:
- The paths in the sealed block (`data/processed/bench/positives_*_TEST.parquet`,
  `regulator_labels_*_TEST.tsv`) match `label_paths()`.
- The seal root is chmod 700 and its sub-directories keep their mode, so reading the files works.
- `.venv/bin/python -m pytest -q tests/test_frozen.py` gives **2 passed**: no drift, and no sealed files are
  back in the tree.
- None of the 37 TRAIN `compare` rows in the tuning log has a NaN CMH bound, so the problem in M2 did not
  affect any tuning decision.

**(5) Resources and crash risk.** The fresh mock took **245 s wall (174 s user) with a peak RSS of 5.96 GB**,
on a 10-CPU machine with 32 GB. TEST is the same size as the mock for HCT116 (1,817 regulators) and slightly
smaller for HEK293T (798 vs 884). Expect about 4 min and under 6 GB. Disk use is negligible: about 50 JSONs
of around 7 kB each.

## Required fixes

### Blocker
None.

### Major (resolve or accept in a decision note before the run; none can be fixed afterwards)

- **M1. Sealed-file pre-check happens after the LOCK.** `check_guards` never checks that `ctx.sealed` and
  the 4 label files exist. If `~/DepMapTR_sealed` were moved, or `DEPMAPTR_SEALED` / `--sealed-dir` were
  set wrongly, `unseal_and_verify` would raise *after* the LOCK commit. The run would be marked FAILED and
  burned, even though no label was read.
  - Fix: before writing the LOCK, add a stat-only check that each label path exists, has an entry in the
    sealed block, and that `ctx.sealed` equals the root recorded in the `# sealed … -> /Users/colmclaude/DepMapTR_sealed`
    line. Optionally also check `p.stat().st_size > 0`.
  - `final_eval.py` is frozen, so this needs a decision note and a re-hash.

- **M2. The paired-bootstrap CMH bound is NaN for many comparisons.** `compare()` computes
  `np.quantile(bc, 0.95)` over replicates that can be ±inf, and one inf − inf gives NaN.
  - In the mock, `boot_cmh_ub90` is NaN for **all 18 HEK293T comparisons** and for HCT116 vs B3, B4 and
    **B5**. The stderr shows "RuntimeWarning: invalid value encountered in subtract".
  - Reproduced on TRAIN fold-2: selected vs B5 has 1 non-finite replicate in 1,000 on HCT116 and 3 on
    HEK293T.
  - "Must beat B5 on pooled **and** CMH OR" is a pre-registered criterion (DESIGN §7). The point estimate
    ΔlogCMH (−1.18 on HCT116) is still reported, but the bound is not.
  - Fix: use `np.nanquantile` over finite replicates and report `n_nonfinite`, or decide now to report only
    the point estimate. This changes frozen `benchmark.py`.

- **M3. The context layer's headline comparisons are not computed on TEST.** `run_all` pairs only
  selected-vs-each. D7's claims about the context layer — vs B5 (Δlog-OR +0.20, not significant) and vs C2
  (24.1 vs 13.9) — therefore get no paired jackknife or bootstrap on TEST. The per-run JSONs do not store
  the per-arm counts needed to compute them afterwards.
  - Fix: add an optional `comparisons: [[RUN_A, RUN_B], ...]` key, e.g. context layer vs B5, vs C2_g2 and
    vs B4e. This changes `final_eval.py` (frozen).
  - If the PI declines, REPORT must state that these are unpaired on TEST.

### Minor (yaml, reporting or hygiene; can be done now or handled in REPORT)

- **m1. NC2 is weakly informative under G1** (see the evidence above). Report a, E[a] and the number of
  regulators that contribute to the top-K next to the OR.
  - Optional, yaml-only: add an NC2 row for the no-G1 ablation:
    `subsets: [..., {name: NC2_noG1, run: a3_r095_s2_w0, subset: nc2}]`.
  - In both NC2 rows, the sign-concordant OR of 3.89 with `n_overlap = 0` is a Haldane artefact. Show NA
    when a = 0.
- **m2. class_iii has no summary row.** DESIGN §6 reports classes (i) and (iii). Class (iii) exists only
  inside each JSON (`metrics.class_iii`) and is missing from summary.md/tsv. Add `class_iii` to `subsets`
  (yaml-only) to get the full metric set.
- **m3. Some reported quantities are in the JSONs but not the summary tables.** Missing from the tables:
  - the sign-concordant Fisher p;
  - the sensitivity label sets (`pos_ueffect`, `pos_p01`);
  - the overlap count a and the number of positives for the pooled OR;
  - the bootstrap CI.

  REPORT must pull these from `results/final/<run>/*.json`. The sign-concordant OR has no
  regulator-permutation p; its significance comes only from Fisher and the within-regulator sign
  permutation. State this in REPORT.
- **m4. The final commit leaves the tuning log uncommitted.** It runs `git add results/final` only, so the
  TEST rows appended to `results/tuning_log.tsv` stay uncommitted and the tree is left dirty (DESIGN §8.3
  says the log is committed). Add `results/tuning_log.tsv` to both commits.
- **m5. Rerun and commit hygiene.**
  - The LOCK guard checks only the working tree and HEAD. Also refuse if
    `git log --all -- results/final/LOCK` is non-empty.
  - Require `git status --porcelain` to be empty (untracked files excepted), so the LOCK commit carries
    nothing else and `git_hash` is not `+dirty`.
- **m6. Score hashes are not recorded before unsealing.** The scores are gitignored, and their sha256s reach
  git only through the per-run JSONs, which are written after TEST is read. Write the `scores.npz` sha256s
  into the LOCK before unsealing. The mock JSONs already provide a committed pre-TEST record, and they
  match the current files.
- **m7. Misleading C1 label.** `c1_mode` reads "score (approximation)" for the selected run. That mode is
  exact for a score with no target factor (BENCHMARK §10). Label it as exact to avoid misreading.
- **m8. Runtime reporting.**
  - `threads` shows `{'cpu_count': 10, 'blas': []}`: BLAS threads are not detected.
  - PC2 has no score wall time.
  - REPORT should join OR and wall time per method in one table ("runtime alongside accuracy").
- **m9. Git identity.** The commit identity is auto-derived from the DHCP hostname
  (`colmclaude@dhcp-…ucd.ie`). Set a repo-local `user.name`/`user.email`, so a hostname change cannot fail
  the LOCK commit. Such a failure would happen before any read, but it would leave an uncommitted LOCK that
  blocks the rerun.
- **m10. selected.yaml hygiene.** When copying DRAFT to `selected.yaml`, update the header comment ("Not used
  by final_eval.py until…"). Optionally relabel the fused run "HCT116 context layer (reported variant; D7)".

## Run checklist

1. Apply or accept M1–M3 with a decision note, then run `record_frozen_change.py` and
   `pytest tests/test_frozen.py`.
2. Copy DRAFT to `results/final/selected.yaml` and commit it.
3. Re-run `final_eval.py --mock --selected results/final/selected.yaml` and check that the output is
   unchanged apart from the intended changes.
4. Make sure the tree is clean and `DEPMAPTR_ROOT` and `DEPMAPTR_SEALED` are unset.
5. Run exactly `.venv/bin/python scripts/final_eval.py` from the repo root. The token check requires
   `argv[0]` to be `scripts/final_eval.py`, so do not use `-m` or a wrapper.
6. Expect about 4 min and under 6 GB.
