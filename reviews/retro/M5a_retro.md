# M5a process retro (fresh reviewer; process only, not science)

## What should we stop doing?
1. Freezing a design before computing the cheap baselines it is judged against. DESIGN s13 (B5 row) says the method "must beat B5 on pooled and CMH OR", yet B5 first appeared at M4 (OR 19.7 vs M1 1.59), after ~1.5 h of review rounds and a frozen, sealed pipeline. Compute DepMap-only baselines at M1/M2.
2. Tuning on an objective already known to reward the wrong thing. D2 tuned pooled OR with a CMH guard. Tier 2 gained +0.79 log-OR while CMH stayed ~1 ("regulator selection, not target ranking"). Pick the objective around the stated scientific claim, or report the greedy path as selection-only.
3. Running multi-round reviews at fixed size regardless of yield. Round 1 = 7 reviewers (~75 KB), round 2 = 5, round 3 = 1 fresh PI. DESIGN.md is 36 KB with 520 lines. Meanwhile code review found a blocker (D5, NaN in 28% of A3 regulators) and D4 found a meaningless threshold. Spend review budget on code and on pre-registered thresholds, and cap design rounds at 2.
4. Writing frozen thresholds without checking the data scale. D4: the |Ueffect-0.5|>=0.1 rule assumed [0,1] but Ueffect is in [-1,1], and it kept 98.6% of positives. A 1-minute distribution check on real data would have caught it before freeze.
5. Letting LESSONS.md go stale. No rule was added after M2 (last commit 57a3621) though M3-M5 produced clear evidence (D3, D4, D5, 22 tuning logs). Log lessons when decisions D3-D5 are made, not at milestones.

## What worked well?
1. Independent narrow-remit reviewers, in parallel. The harness review found no blockers but flagged a "run once" crash hole and tuning-log bypass; the model review found the A3 NaN blocker. Both were fixed before the seal (commits 1ed1efa, dfe403c).
2. Freeze-then-seal discipline. The TEST labels were moved out of the tree with hashes committed (57f2de3), and the split is by arm, set before tuning.
3. Pre-registered, greedy tiers with a strict adoption rule (fold consistency + jackknife lb90). It rejected 5 of 7 tiers, including the tempting focal-arm gain (+0.68 but fold1 -0.09 / fold2 +1.52).
4. Honest negative reporting. RESEARCH_TREE states that the incumbent is far below B5 and fails CMH. Hub baselines and permutation controls (B0, NC1, B6) exposed the hub dominance.
5. Parallel data/harness/model agents with docs/INTERFACES.md fixed up front. M3 took about 70 minutes with no interface rework visible in the log.

## Proposed new LESSONS (one line each)
- Compute the trivial DepMap-only baselines (hub, expression, co-expression) before freezing the design; the freeze criteria should cite their numbers. Evidence: B5 19.7 vs M1 1.59 appeared only at M4.
- Print the distribution of every real-data column a frozen threshold uses, before the freeze. Evidence: the D4 Ueffect rule kept 98.6% of positives.
- Add a LESSONS line at each D-number decision that cost effort. Evidence: none added since M2 despite D3-D5.

## Prune candidates
- #1 (check environment) and #2 (read README/file sizes) are generic one-off M0 steps; merge into one "pre-flight" line or drop.
- #5 (continue reviewers via SendMessage) is an optimisation with no measured evidence; keep only the "fresh reviewer last round" half.
- #4 and #6 are valuable but science-specific; keep them, but tag them as data-audit rules.
