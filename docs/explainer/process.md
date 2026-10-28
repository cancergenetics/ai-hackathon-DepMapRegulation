# How the research was actually done

**Sources.**
- `git log --date=iso` (UTC).
- RESEARCH_TREE.md and LESSONS.md.
- reviews/retro/M5a_retro.md and reviews/retro/M7_retro.md.
- The main session transcript. It was searched only for the timestamps of human messages and agent launches; it
  was not read in full.

## 1. Numbers at a glance

| item | count / span | source |
|---|---|---|
| Wall-clock, brief to project close | 17:58 → 23:28 on 2026-10-27, **≈ 5 h 30 min** | brief pasted 17:58:15; first commit `e31c6f7` 18:03; close `c27a904` 23:28 |
| Commits | 40 during the project, plus 1 the next morning (README/disclaimer, `a0f02e0`) | git log |
| Subagents launched during the project | **19**: 14 Opus, 5 Sonnet | main transcript |
| Follow-up messages to running agents | 16 | main transcript |
| Design-review runs | 13: 7 (round 1) + 5 (round 2, the same agents continued) + 1 (round 3, fresh) | reviews/round1–3 |
| Code / protocol reviews | 2 component reviews (harness, model) + 1 pre-final review | reviews/code |
| Other agents | 1 literature, 1 data, 1 harness builder, 1 model builder, 1 ideas scientist, 1 report verifier, 2 retros | |
| TRAIN tuning configurations | **29 of a 30 cap**: 23 pre-registered (M1 + tiers 1–7) + 6 post hoc (D6) | results/tuning_log.tsv |
| TRAIN sensitivity / ablation runs (reporting only) | 7 | RESEARCH_TREE D8 |
| TEST runs | **1** LOCKed run: 55 evaluations of 23 score runs over 2 lines, 4.5 min | results/final/ |
| Tests | 187 (SIM known-answer, OLS agreement, metrics) | REPORT §8 |
| Review text produced | ≈ 75 KB (round 1), 28 KB (round 2), 6 KB (round 3), 41 KB (code reviews), 19 KB (ideas) | wc -c |

Three more agents were launched on 2026-10-28 for this write-up, including the one that wrote this page. They are
not counted above.

## 2. Phases (wall-clock, from commits)

| phase | span | duration | what happened |
|---|---|---|---|
| M0 setup | 17:58–18:03 | 5 min | PROMPT.md saved verbatim, venv, research tree |
| M1 literature + data | 18:02–18:33 | 31 min | literature and data agents ran in parallel |
| M2 design + 3 review rounds | 18:33–18:55 | **22 min** | PI wrote v1; 7 reviewers; v2; 5 reviewers; v3; 1 fresh reviewer; v3.1 frozen |
| M3 build + code review + seal | 18:56–20:05 | 69 min | data, harness and model agents in parallel against docs/INTERFACES.md; two independent reviews; D5 fix; TEST sealed |
| M4 baselines | 20:05–20:08 | 3 min | M1 1.59 vs B5 19.7 (TP5) |
| M5a pre-registered tiers | 20:08–20:47 | 39 min | G1 and A3 adopted (TP6) |
| M5 post-hoc tier (D6/D7) | 20:47–22:02 | **75 min** | retro, ideas scientist, I1–I6, C1/C2; nothing added to the CN-only selection; context layer split off (TP7) |
| M6 preparation | 22:02–23:08 | 66 min | D8 sensitivity analyses, mock final run, pre-final review, D8a fixes |
| M6 TEST run | 23:08–23:12 | 4.5 min | selection commit, then LOCK commit, then one run |
| M7 report | 23:12–23:28 | 16 min | REPORT draft, independent verification, corrections, final retro |

## 3. Who did what

- **PI (main agent).** The PI did the following itself:
  - wrote DESIGN.md v1 and every revision, including all 3 × §13 disposition tables;
  - wrote all decision notes D0–D8a and the research tree;
  - wrote LESSONS.md, the first REPORT draft and the corrections;
  - built and ran the tuning loop (`scripts/tune_tier.sh`, the compare runs);
  - made every adopt/reject call;
  - ran the seal, the mock run and the final run;
  - made every commit, from its own shell. Two of them (the LOCK commits `603a7f5` and `e9085c7`) were written by
    `final_eval.py`, which the PI ran.

  Its own tool use was mostly shell: about 140 Bash calls, 35 edits and 10 file writes over the session.
- **Literature agent (Opus).** Wrote LITERATURE.md with 59 references, verified on the web, and marked unverified
  numbers. It later applied the round-1 corrections (T5): ERBB2 is on 17q, the CausalBench conclusion was updated,
  and a Warmerdam caveat was added.
- **Data agent (Opus).** Wrote DATA.md, DATA_PROFILE.md and the data scripts and tests. It inferred the release as
  26Q3. Its data profile produced D1 (TP2). It later ran T1–T4: TP53, RNA technical covariates, the S2 panel and the
  frozen flags.
- **Design reviewers.**
  - Round 1 had 7 agents, all independent and in parallel: external PI, statistics, data, domain biology (all Opus);
    prior work, algorithms/compute, AI-agent failure modes (Sonnet).
  - Round 2 continued 5 of those agents on v2.
  - Round 3 was one *fresh* external PI on v3.
  - Most of the evaluation design (TP3) and the specification of the gate, attribution and sign (TP4) came from
    these reviewers.
- **Harness builder and model builder (Opus).** Built in parallel against fixed interfaces. They received follow-up
  work for D3, D5, D6, D7/D8 and D8a.
- **Code reviewers (2, Opus).** Each re-implemented the component independently.
  - Harness review: no blockers, 2 procedural majors (a hole in the "run once" guard; tuning-log bypass).
  - Model review: core model reproduced to 5e-8. It found the A3 blocker: NaN scores for 28% of pairs with
    w_co > 0 (D5).
- **Ideas scientist (Opus, fresh).** Diagnosed the CN signal on TRAIN only:
  - within-regulator AUROC 0.507;
  - about 94% of the log-enrichment comes from regulator allocation;
  - t-rows track the arm signature at r = 0.87.

  It proposed I1–I6 and the C1/C2 controls (TP7). This is the most consequential single agent contribution after
  round-1 review.
- **Pre-final reviewer (Opus).** Ran an independent mock of the final run from a clean clone. Gave a conditional GO
  with 3 fixes (TP8).
- **Report verifier (Opus).** Checked about 120 numbers. Found 4 errors and 3 over-claims.
- **Retros (2, Sonnet, fresh).** Answered "stop doing / worked well" at M5a and M7. LESSONS.md was pruned from both.
- **Saved by the PI.** The harness blocked subagents from writing report files, so the PI wrote the files
  `reviews/report/report_verification.md` and `reviews/retro/M7_retro.md` as condensed versions of what the agents
  replied. They are not verbatim agent output.

## 4. Where exploration was productive

- **Parallel, narrow-remit reviewers.** At least 4 of 7 round-1 reviewers independently flagged pair dependence. Four
  independently asked for an arm-level split. Issues raised by several reviewers turned out to be real (LESSONS #4).
- **Fixed interfaces before parallel builds.** Data, harness and model were built in about 70 minutes with no
  integration rework (LESSONS #6).
- **Independent re-implementation in code review.** It caught three problems before any evaluation:
  - the A3 NaN blocker;
  - the LOCK crash hole;
  - the tuning-log bypass.

  No design reviewer could have seen these (LESSONS #3, #7).
- **The strict adoption rule.**
  - It rejected the unstable focal (+0.68) and I4 (+0.21) gains.
  - Neither looked better on TEST: focal 4.50 vs 3.85, n.s.
  - The G1 gain it did adopt replicated on TEST: +0.97, CI 0.30–1.64.
- **The ideas scientist and the C1/C2 controls.** They turned "we are far below B5" into a precise statement: CN
  helps choose regulators and directions, not targets. They also stopped a +1.57 pooled gain from being reported as
  CN knowledge.
- **The one-shot TEST protocol.** It followed the sequence seal (`57f2de3`), selection commit (`62e67b9`), LOCK
  (`603a7f5`), run (`e9085c7`), with a mock dry run first. The pre-final review found 3 real fixes in about 13 min
  (LESSONS #8).

## 5. Where it went in circles or wasted effort

- **Design review was heavy in tokens, light in time.**
  - The 13 review runs produced about 110 KB of text. DESIGN.md grew to 36 KB / 520 lines.
  - Round 3 had to fix contradictions that round 2 introduced: the sign rule differed between §1 and §6, c₀ was
    undefined, and round-2 items were silently dropped.
  - The M5a retro called this "~1.5 h of review rounds". The commits show 22 min (18:33 → 18:55); the 1.5 h figure
    fits M1→M4 as a whole.
  - The retro's main point stands: the blockers came from code review, not design review (LESSONS #3).
- **Baselines were specified but not computed before the freeze.** B5 was added in round 1, and the frozen design
  says the method "must beat B5". Its value (19.7, against M1's 1.59) appeared only at M4, after the seal. The
  method never met that requirement (TRAIN Δlog-OR −1.37) (LESSONS #1).
- **The objective was chosen literally from the brief.**
  - Statistics #11 proposed tuning on CMH. The PI kept pooled OR because the brief asks for enrichment.
  - As a result, every adopted gain was regulator allocation (G1 ΔCMH −0.05).
  - The one variant that improved within-regulator CMH (I3, +0.42) was rejected.
  - Both retros list this as a "stop doing" (M5a #2, M7 #2).
- **A proposal was dismissed, then returned post hoc.**
  - Round-1 prior-work reviewer #5 proposed a DepMap-only target-expression prior and predicted it would be the
    largest OR gain.
  - DESIGN §13 marked it "covered by B3/G1". That is inaccurate: B3 and G1 cover only the regulator side.
  - The idea returned about 2.5 h later as the post-hoc I1 fusion.
  - Pre-registering it would have avoided the post-hoc label on the context layer.
- **The D6 tier took 75 minutes and 6 configurations and changed nothing in the selected method.** It did produce the
  S5 layer and the decisive diagnosis, so this is partly "wasted" and partly the most informative hour of the
  project.
- **Several mis-specifications were caught late:**
  - D4: the Ueffect threshold was written for the wrong scale and kept 98.6% of positives.
  - D6a: C1 silently ran in a biased approximate mode.
  - I2's adoption rule could not be satisfied by construction.
  - The NaN CMH bound.
  - The missing sealed-file preflight.

  None of these had a known-answer test before first use (LESSONS #11).
- **Record-keeping lagged.**
  - LESSONS.md got no new lines between M2 and the M5a retro, nor for D6–D8a.
  - REPORT.md had 4 errors and 3 over-claims (Fisher p 2×10⁻³⁶ shown next to the pre-registered p 0.068), and the
    research tree repeated them until the M7 retro (LESSONS #12).
  - Some residual over-claims remain in the tree today; see insights.md, "Inconsistencies".

## 6. Where human input changed direction

I checked the main session transcript for every message typed by the human. The only human inputs were:

1. **2026-10-27 17:55.** Session setup: `/effort high` and `ls`.
2. **2026-10-27 17:58.** The research brief (PROMPT.md, verbatim). The brief carried substantive priors that shaped
   the method:
   - each gene's own cis CN must be modelled;
   - trans means different chromosomes;
   - correlated arm events can confound;
   - a suggested covariate list;
   - the evaluation universe;
   - "CN loss should correspond to knockdown";
   - hold out perturbed genes before tuning;
   - "maximise enrichment (odds ratio and significance)".

   The last point drove the pooled-OR objective (see §5).
3. **2026-10-28 11:41–11:45, after the project closed.** "Commit the code to a new repository", then (after the PI
   asked about visibility) "private repository, do not commit any of the data files, add the following disclaimer".
   The result was `a0f02e0` (README with the AI-hackathon disclaimer; npz/csv ignored).
4. **2026-10-28 12:00.** The request for this write-up. The human answered one question ("Both") on the deliverable.

**What this means:**
- Between 17:58 and 23:28 there were **no human messages**. All research decisions (D0–D8a, every adoption call,
  the stop decision) were made by the PI agent.
- The commit history agrees: all 40 project commits came from the PI agent's shell. The author name "Colm Claude"
  is the machine's git identity, not a sign of human edits.
- No decision note cites human input.
- The two `AskUserQuestion` calls both came after the project closed (repository visibility; the deliverable).

## 7. Why it stopped

- The PI stopped at 23:28. Its reason: TEST is LOCKed, so further method changes cannot be validated without a new
  hold-out.
- One of the 30 configurations was left unused.
- Ideas remained:
  - near-misses that need a fresh hold-out: I3, focal, I4;
  - ideas that need new data: focal-event instruments, proteomics, joint models.

  The stop was therefore a protocol boundary, not exhaustion of ideas.
- The brief asked to continue "until you genuinely run out of promising ideas or close to the time limit". Arguably,
  more budget could have gone into label-free diagnostics before the single TEST run, but the PI chose not to spend it.
