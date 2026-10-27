# LESSONS — process rules for running this project efficiently

Each rule is one line with its evidence. Pruned at each milestone (last: M5a, after the fresh-agent retro
reviews/retro/M5a_retro.md).

## Planning and design
1. Compute trivial label-free baselines (hub, detectability, co-expression) before freezing the design. Evidence: B5 (DepMap-only hub prior) reached OR 19.7 vs M1 1.59, but only surfaced at M4, after three review rounds.
2. Print the real-data distribution behind every frozen threshold before freezing. Evidence: D4. The Ueffect rule assumed a [0,1] scale and kept 98.6% of positives, so it was meaningless.
3. Cap design review at 2 full rounds plus 1 fresh-eyes check; spend the saved effort on code review. Evidence: rounds 1–2 produced ~75 KB of reviews, while code review found a blocker (D5 NaN propagation) that no design reviewer could have seen.
4. Run independent reviewers in parallel, each with a narrow remit; issues raised by several of them are almost certainly real. Evidence: 4 of 7 round-1 reviewers flagged pair dependence in Fisher p; 4 flagged 5 Mb split leakage.

## Data audits (science-specific)
5. A second agent audits the data agent's numbers before the design relies on them. Evidence: S1 silently excluded HCT116 itself (M2 audit E2).
6. Sex chromosomes break naive CN statistics (male X ≈ log2 −1): sex-centre X before thresholding. Evidence: 95% of X genes passed the raw variable-CN rule.

## Execution
7. Fix interfaces (docs/INTERFACES.md) before spawning parallel builders. Evidence: the data, harness and model agents built M3 in parallel in ~70 min with no integration rework.
8. Independent code review with the reviewer's own re-implementation before any evaluation. Evidence: the harness review caught a crash hole in final_eval LOCK; the model review caught 28% NaN scores (D5). Both were fixed before the seal.
9. Physically seal TEST labels (outside the tree, unreadable, hashed) rather than relying on flags. Evidence: the round-2 AI-failure review showed that flag-based guards were bypassable.
10. A strict paired adoption rule (both folds > 0 and paired lower bound > 0) prevents chasing noise in hub-dominated benchmarks. Evidence: focal arm handling Δlog-OR +0.68 came from one fold (−0.09 / +1.52) and was rejected.
11. When the tuning objective can be gamed by regulator selection, track a within-regulator metric (CMH, degree-stratified permutation) next to it at every tier. Evidence: G1 gained +0.79 log-OR with CMH unchanged (≈1).
12. Add a LESSONS line at every decision note (D-number) that cost or saved effort. Evidence: no lessons were added between M2 and M5a despite D3–D5.
