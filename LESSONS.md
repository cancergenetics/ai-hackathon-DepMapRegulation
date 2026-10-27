# LESSONS — process rules for running this project efficiently

Each rule is one line with its evidence. Pruned at each milestone. Last pruned at M7 (project close), after
the fresh-agent retros reviews/retro/M5a_retro.md and reviews/retro/M7_retro.md.

1. Compute label-free baselines (hub, detectability, co-expression) before freezing the design, and cite their numbers in the freeze criteria. Evidence: B5 (DepMap-only hub prior) at OR 19.7 vs M1 1.59 surfaced only at M4.
2. Print the real-data distribution behind every frozen threshold before freezing. Evidence: D4. The Ueffect rule assumed a [0,1] scale and kept 98.6% of positives.
3. Cap design review at 2 rounds plus 1 fresh check, and spend the saved effort on code review. Evidence: code review found 28% NaN scores in A3 (D5), which no design reviewer could have seen.
4. Run reviewers in parallel with narrow remits; issues raised by several of them are real. Evidence: 4 of 7 round-1 reviewers independently flagged pair dependence of Fisher p.
5. A second agent audits the data agent's numbers before the design uses them; sex-centre X in CN statistics. Evidence: S1 silently dropped HCT116, the benchmark's own line (audit E2); 95% of X genes passed the raw variable-CN rule.
6. Fix interfaces (docs/INTERFACES.md) before spawning parallel builders. Evidence: data, harness and model were built in parallel in ~70 min with no integration rework.
7. Code review uses the reviewer's own re-implementation, both before the first evaluation and before the final run. Evidence: these reviews caught the LOCK crash hole, the NaN propagation, and 3 pre-final fixes (sealed-file preflight, NaN CMH bound, missing comparisons).
8. Seal TEST physically, commit the selection, commit the LOCK, then do one run, after a mock dry run on the same code path. Evidence: flag-based guards were bypassable; the protocol gave one clean run (62e67b9 → 603a7f5 → e9085c7).
9. Use a strict paired adoption rule (both folds > 0 and paired lower bound > 0) on hub-dominated benchmarks. Evidence: focal arm handling (+0.68) came from one fold (−0.09 / +1.52) and was rejected; on TEST it was not significantly better.
10. Pair every pooled objective with within-regulator controls (CMH, pair-shuffle C1, permuted-CN fusion C2) at every tier, and run them before adding a new tier. Evidence: G1 gained +0.79 with CMH ≈ 1; the I1 fusion's +1.57 vanished under C1 (D7).
11. Every new metric or control gets a toy-data known-answer test before its first real use; never write an adoption rule that a variant cannot satisfy by construction. Evidence: D6a (C1 silently in approximate mode); the late NaN CMH bound; the D6 rule was unsatisfiable for I2.
12. Draft claims from script-generated tables, show only the pre-registered p next to headline numbers, and have a fresh agent verify every summary document (report and research tree). Add a LESSONS line in the same commit as each decision note. Evidence: REPORT.md had 4 errors and 3 over-claims (p 0.068 next to Fisher 2e-36); the tree repeated them until the M7 retro.
