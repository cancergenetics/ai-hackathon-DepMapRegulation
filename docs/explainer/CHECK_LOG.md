# Independent check of the explainer materials (2026-10-28)

Scope: README.md (top part; Disclaimer untouched, verified byte-identical), docs/site/{METHOD,JOURNEY,INSIGHTS}.md,
docs/site/report.html (+ rebuilt index.html), docs/briefing/BRIEFING.md, docs/explainer/{turning_points,insights,
process,dead_ends,worked_example}.md, journey_timeline.json. Ground truth: results/final/summary.{tsv,md},
results/final/*/hct116_TEST.json, FACTS.md (incl. Correction D9), RESEARCH_TREE.md, LITERATURE.md,
reviews/ideas/M5a_ideas.md, docs/DATA_PROFILE.md, worked_example.json, git log.

## What was checked and found correct

- Headline and paired numbers (3.85, CI 1.33–11.2, p 0.068; G1 +0.97 CI 0.30–1.64; NC1 +1.10 CI −0.66–2.85; B1 +1.22
  CI −0.48–2.92; C2 vs selected −1.01 CI −1.99 to −0.04; context vs B5 −0.22, vs C2 +0.55 CI 0.07–1.03; sign 0.718 /
  0.647 / p 0.012; C2 sign 0.585 / 0.532 / c₀ 0.514 / p 0.001 / n 354) re-read from summary.md and the TEST JSONs.
- TRAIN tuning path, diagnostic numbers (0.87, 91 %, 1,308, 1,701, 0.0006→0.0124, 3.5 %/1.7 %, AUROC 0.70/0.511/0.498,
  log ABF 0.385 z², 28.5 %, 98.6 %, 47 %, 1,255) against M5a_ideas.md, DATA_PROFILE, BENCHMARK, code review.
- Worked-example numbers against worked_example.json.
- Literature credits against LITERATURE.md keys [F1 K1 D2 T1 S4 D3 S5 T2 G2 B2 W1 W4 A1 L2 F2] plus Replogle 7.5 %,
  Huang 75.4 %/51.5 %, Schukken (flagged unverified, as LITERATURE.md does), Warmerdam OR 12.4 on 41 pairs,
  "ORs about 1–2", Goh siRNA validation. Wakefield's Bayes factor is correctly labelled as a reviewer-agent idea not in
  LITERATURE.md; LOCO correctly labelled "GWAS practice, no citation" (LITERATURE.md itself recommends it).

## Change log

### Overclaiming / underclaiming
1. README "In short": "copy number tells you which regulators…" → "the results suggest that copy number says something
   about which regulators… (gate replicated, overall enrichment not significant); no information about which targets
   was detected".
2. README results table: added "not significantly above a shuffled-copy-number control (+1.10, CI −0.66 to 2.85)".
3. README steps: step 1 notes expression-PC harm is literature-only; step 3 notes attribution not tested on held-out
   data; step 5 adds "held-out OR 18.4, no better than a label-free hub prior (22.8)".
4. METHOD.md + report.html, "method in three sentences": sentence 3 was a result, not a method step, and said the
   method "identifies which regulators". Now sentence 3 = score/rank/direction; a separate "What it showed"
   paragraph gives p 0.068, n.s. vs shuffled CN, the gate as the one replicated result, "suggests … which
   regulators", no detectable target information, direction not established.
5. METHOD.md + report.html Results lead: "all of it comes from choosing regulators" → "what enrichment there is appears
   to come from choosing regulators".
6. METHOD.md + report.html Step 5 / context layer: "+0.55 vs C2 … so real copy-number evidence helps choose regulators"
   → "which suggests …; post hoc, not multiplicity-corrected".
7. report.html pipeline "Output / Use it for": added "(weakly supported: gate's gain replicated, overall enrichment
   not significant)".
8. "The signal is which regulators" figure (report.html inline SVG, briefing PNG/SVG via make_figures.py, figures.json,
   BRIEFING Fig. 6, METHOD alt text): title → "What signal there is: which regulators…"; panel (a) "Remove the CN
   evidence: OR falls" → "Permute CN: OR falls (CN-only: n.s.)"; footnote now gives paired Δlog-OR (CN-only +1.10,
   CI −0.66 to 2.85, n.s.; context layer +0.55, CI 0.07–1.03); captions say the CN-only drop is not significant.
9. BRIEFING §8.1: "Four independent lines … Permuting copy number does hurt" → "Four lines … the last not
   significant on its own"; NC1 bullet now states the paired difference is n.s.
10. BRIEFING §7 finding 4, §4.3, §10: "+0.55 … adds something" → "suggests … (post hoc)"; attribution "never ablated
    on TEST"; "outscore every copy-number ranking" → "outscore the copy-number-only ranking (3.85)" (the context layer
    beats C2); "What was not shown" now lists enrichment over NC1, advantage of context layer over hub prior, and
    held-out value of attribution.
11. INSIGHTS.md + report.html Insights tab: main lesson and "Which regulators" highlight softened ("if they say anything
    here … weakly supported"; "more likely to have", with CI and the n.s. headline); "no pair-level knowledge at all"
    → "no detectable pair-level knowledge".
12. JOURNEY.md + report.html: TP7 title "no target-level information" → "no detectable …"; TP8 outcome "copy number
    identifies which regulators" → "may help identify … (resting mainly on the gate's replicated gain)"; why-section
    "a diagnosis showed no target-level information" → "found no detectable …".
13. Sign-concordance figure footnote and report.html caption stated the "correlated modules" explanation as fact →
    "Untested explanation: …" / "A likely but untested reason: …" (make_figures.py, SVG/PNG, report.html).

### D9 consistency (direction not established)
14. docs/explainer/insights.md #6 heading "CN-predicted direction agrees modestly…" → "… is not established
    (correction D9)"; C2 now quoted against its within-regulator null (0.532), not only c₀; "REPORT.md does not
    discuss this" → corrected under D9.
15. turning_points.md TP8 outcome "and modestly their direction" → withdrawn under D9; "replicates TRAIN +0.79" →
    "+0.91 for the same removal contrast; +0.79 when added to M1" (FACTS discrepancy 3).
16. journey_timeline.json TP8 one_line "modest direction" → withdrawn (D9). JSON re-validated.
17. worked_example.md "Direction is right more often than not … p = 0.012" → "whether that beats chance is not
    established (C2 passes the same test, p = 0.001)".
18. process.md "CN helps choose regulators and directions" → "may help choose regulators … direction claim withdrawn
    (D9)"; "Neither [focal nor I4] looked better on TEST" → I4 was never run on TEST.
19. JOURNEY.md / report.html "corrections" appendix updated to say those notes were corrected, and that the
    RESEARCH_TREE.md status summary still carries the old wording.

### Numbers / precision
20. METHOD.md + report.html "Cost: peak footprint 6–8 GB" → "6.4–8.4 GB memory footprint (macOS max RSS 14 GB)"
    (FACTS §7; 8.4 > 8 and max RSS was omitted).
21. Worked example: PIK3C2B "α ≈ 0.002" is the median over targets (0.0002 for MYBL2) and MDM4's "14 %" is a median
    (6.6 % for MYBL2) — labelled in report.html and METHOD.md.

### Jargon
22. METHOD.md + report.html terms table: added Δlog-OR and per-regulator AUROC / CMH (used undefined in results).
23. Pseudocode intro (METHOD.md, report.html) now defines the comment ids G1, A3, NC1, I1.
24. JOURNEY.md "LOCO kept" → "the leave-own-chromosome-out arm adjustment kept"; INSIGHTS "(D9)" → "(decision note
    D9)"; README "AUROC 0.500 (0.5 = chance)".

### HTML hygiene
25. Removed the three Google Fonts `<link>` tags (the only external fetch); CSS font stacks fall back to system fonts
    (Arial Narrow/Helvetica, Georgia, Menlo). Re-add if the page is published as an Artifact and the fonts are wanted.
26. Rebuilt docs/site/index.html with `.venv/bin/python scripts/figures/build_site.py`. Checks on report.html: no
    doctype/html/head/body tags; 290 `href="#…"` anchors, 0 unresolved, 0 duplicate ids; no external src/href/url();
    55 colour/layout tokens on bare `:root`, dark overrides under `@media (prefers-color-scheme: dark)
    :root:not([data-theme="light"])` and `:root[data-theme="dark"]`; all `var(--…)` defined; inline JSON parses; inline
    JS (19.3 kB) passes a JavaScriptCore syntax check (`osascript -l JavaScript`, `new Function`; node not installed;
    a deliberately broken file fails as expected). Headless Chrome render of the Method tab inspected: three sentences,
    "What it showed", glance row and terms table render correctly with fallback fonts.
27. Regenerated docs/briefing/figures with make_figures.py: only sign_concordance.{svg,png}, what_the_signal_is.{svg,png}
    and figures.json changed (all other SVGs byte-identical).

## User checklist (after edits)
- (a) First screen explains the method in three sentences: **yes** — README (one sentence + five steps) and Method tab
  (one sentence + three method sentences, results separated).
- (b) Every step has a "why": **yes** (Steps 1–5 and Evaluation each carry Why / Without it).
- (c) Claims tied to evidence, unverified ones marked: **yes**, after edits 1–13 (literature-only, post-hoc,
  untested-explanation and not-ablated items now flagged).
- (d) Journey explains why the method looks as it does: **yes** — "Why the method looks the way it does" maps each
  feature to a turning point; the context-layer reason was clarified (edit 12).

## Open issues (not fixed)
- RESEARCH_TREE.md (ground truth, not edited) still says "a modest excess of correct effect directions (p 0.012)" in
  its status summary (line 10), "(ii) the direction of effects" (line 30) and "[!] M1 has direction information"
  (M4). Should be corrected by the PI to match D9.
- REPORT.md was not re-checked here beyond what FACTS.md records.
- Figure-panel internal labels in the journey timeline (raw commit messages: D6a, I1–I6, G1+A3…) and the tuning-path
  footnote (M1, G1, A3, C1) are not expanded inside the figures; the surrounding captions/text define them.
- The sign-concordance figure subtitle still says "(D9)" without expansion (caption explains it).
- docs/explainer/dead_ends.md quotes jackknife lower bounds where FACTS uses bootstrap ones (already disclosed in
  JOURNEY.md); not changed.
