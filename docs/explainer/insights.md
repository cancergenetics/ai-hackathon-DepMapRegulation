# What was learned

- Each insight is ranked within its category by generality first, then surprise.
- Format: one sentence; **evidence** (number and source file); **arose at** (turning point in
  [turning_points.md](turning_points.md)).
- **[not shown]** marks something believed but not demonstrated in this project.
- Unless stated, numbers are HCT116 OR at the top 0.1% of pairs. TRAIN numbers come from `results/tuning/` or
  `results/m4/`; TEST numbers from `results/final/summary.md`.

---

## Top three in each category

**Scientific**
1. **Which regulators, not which targets.** Cross-line copy-number associations show *which* regulators have trans
   effects that CRISPRi can detect: those whose own mRNA strongly tracks their CN. They show nothing detectable about
   *which* of a regulator's targets respond.
2. **A DepMap trans row is mostly its chromosome arm's signature**, not something specific to the gene.
3. **In a single-line Perturb-seq screen, regulator hubness and target expression depth largely decide which pairs
   count as hits.** A label-free DepMap prior outscores every CN-based ranking.

**Methodological**
1. **Split pooled enrichment into its parts.** On a hub-dominated perturbation benchmark, separate regulator
   allocation, target detectability and pair-specific evidence, using label-free hub priors and within-regulator
   controls. A large pooled gain can contain zero pair knowledge.
2. **Tune on a metric that measures the claim.** Greedy tuning on pooled OR adopts only allocation gains and rejects
   pair-specific ones.
3. **Give every control and metric a known-answer test before its first real use.** A bugged control nearly
   produced a false discovery.

---

## 1. Scientific and domain insights

1. **DepMap CN carries no detectable information about which targets of a given regulator respond to its
   knockdown.**
   - Evidence: within-regulator AUROC 0.507 on TRAIN (333 regulators; reviews/ideas/M5a_ideas.md §1.1) and 0.500 on
     TEST (179 regulators). Regulator-stratified CMH 0.94 on TEST (CI 0.77–1.15, which bounds any lift at about
     15%). Pair-shuffle control C1: 3.76 vs 3.85.
   - Focal (within-arm) CN does no better: AUROC 0.511.
   - Arose at: TP7, confirmed at TP8.
2. **A regulator's cis CN→mRNA instrument strength predicts whether its knockdown has detectable trans
   consequences.** This is the Mendelian-randomisation logic, and it works even though the instrument is measured
   across all lineages rather than in the screened line.
   - Evidence: the G1 gate (F ≥ 100) gives TRAIN Δlog-OR +0.79 (folds 0.65 / 1.12), monotone in τ. On TEST: 1.46 →
     3.85, Δ +0.97 (CI 0.30–1.64).
   - The gain is allocation, not ranking: ΔCMH −0.05.
   - Arose at: TP4 (design) and TP6 (adoption); replicated at TP8.
3. **A regulator's LOCO trans-association row is mostly the signature of its chromosome arm.**
   - Evidence: median r = 0.87 between a G1 regulator's t-row and its arm-only row; 91% of rows exceed 0.7.
   - Redundancy: 54 regulators on 1q share 1,308 top-K targets.
   - The arm-only score itself does not enrich: B6 TRAIN 0.97, TEST 0.61 (M5a_ideas §1.3).
   - Arose at: TP7.
4. **In single-line CRISPRi Perturb-seq, which pairs count as hits is decided mainly by regulator hubness and target
   expression depth.**
   - Label-free DepMap hub prior B5: OR 19.7 (TRAIN) and 22.8 (TEST), against 5.0 and 3.85 for the CN tool.
   - Within a regulator, target expression in HCT116 alone reaches AUROC 0.70.
   - Even the context fusion built on *permuted* CN (C2, 10.6) beats the real CN-only tool on TEST: Δ −1.01
     (CI −1.99 to −0.04).
   - Arose at: TP2 and TP5; quantified at TP7.
5. **Whether a regulator is expressed in the screened line sets its hit yield, by about 20×.**
   - Evidence: among G1 regulators, the per-pair hit rate rises from 0.0006 to 0.0124 across HCT116 expression bins
     of A (M5a_ideas §1.3).
   - Non-expressed HCT116 perturbations average 5.5 trans hits, against 51.0 for expressed ones
     (DATA_PROFILE §10).
   - A *population* instrument is not a *context* instrument. MYCN passes G1 through neuroblastoma amplicons, yet is
     unexpressed in HCT116; it took 1,701 top-K slots for 1 hit.
   - Arose at: TP2, TP7.
6. **CN-predicted direction agrees modestly with knockdown direction. The evidence is weak.** **[not shown
   robustly]**
   - Evidence: TEST concordance 0.718 against a within-regulator null of 0.647, p = 0.012, n = 131. Only 6
     regulators have ≥ 5 overlaps.
   - Not significant in the S1 sample set (p 0.13) or in the no-split-sibling subset (p 0.48).
   - The C2 negative control also "passes" this test (0.585 vs 0.514, p = 0.001); see Methodological #8.
   - Arose at: TP8.
7. **LOCO t-statistics have a negative skew that grows with target expression. Part of it is a TPM compositional
   artefact; the rest is an own-chromosome aneuploidy response.**
   - Skew: t ≤ −3 is twice as common as t ≥ 3 (3.5% vs 1.7%).
   - Mean t falls from −0.075 to −0.311 across target-expression quintiles.
   - Renormalising TPM to exclude chr(A) moves the chr1 row mean from −0.12 to +0.09. It does not change the ranking
     (I5 Δ +0.05). Source: M5a_ideas §1.5.
   - Arose at: TP7.
8. **The X-Atlas/Orion HEK293T trans hits are dominated by perturbations whose own gene is not expressed. They are
   probably artefacts.**
   - Evidence: 83.0% of HEK293T trans hits, against 12.7% in HCT116; their up/down balance also differs (59% up)
     (DATA_PROFILE §10–11).
   - The CN tool fails on HEK293T TEST: OR 0.49, 6 overlaps.
   - Arose at: TP2.
9. **The CN signal is not TF-centred.** **[not shown: low power]**
   - Evidence: RBP regulators OR 3.56 (30 top-K hits); TF regulators 0 hits.
   - The displayed TF and kinase ORs (6.7, 42.9) are Haldane artefacts of empty cells (REPORT §4).
   - Arose at: TP8.

## 2. Methodological insights that generalise

1. **Split any pooled enrichment into three parts — regulator allocation, target detectability and pair-specific
   evidence — before crediting a method.**
   - Evidence: the I1 fusion gained +1.57 log-OR and passed the strict adoption rule, but the exact within-regulator
     shuffle C1 gave the same OR (24.5 vs 24.1; Δ −0.02, CI −0.09 to 0.06).
   - About 94% of the incumbent's log-enrichment was regulator allocation (M5a_ideas §1.2).
   - Arose at: TP5, TP7.
2. **Choose the tuning objective around the scientific claim, not the headline metric.**
   - Evidence: tuning on pooled OR adopted only allocation gains (G1 ΔCMH −0.05). It rejected the only variant that
     improved within-regulator CMH (I3 shrunk-IV, ΔCMH +0.42, pooled lb90 −0.10).
   - Statistics reviewer #11 had proposed a CMH objective at round 1.
   - Both retros list the choice as a mistake.
   - Arose at: TP3, TP6, TP7.
3. **Every control, metric and threshold needs a toy known-answer test before first real use.**
   - Evidence: an approximate C1 (a loader bug, D6a) reported a "pair-specific gain" of +0.063 (lb90 +0.002); the
     exact C1 gave −0.017.
   - The frozen Ueffect rule assumed the wrong scale and kept 98.6% of positives (D4).
   - I2's adoption rule could not be satisfied by construction.
   - Arose at: TP7 (D4 at TP3).
4. **Stratified "pair-specific" statistics must be stratified on the context that drives detectability.**
   - Evidence: CMH stratified by pan-DepMap expression deciles still gave B4e, a target-only score with no pair
     information, CMH 3.8 (TRAIN) and 4.1 (TEST).
   - HCT116-stratified deciles reduce but do not remove this (3.4 TRAIN).
   - Arose at: TP5, TP7.
5. **Hold out CN-derived predictions by chromosome arm, and budget for the precision this costs.**
   - Gain: a TEST regulator's maximum CN r with any TRAIN regulator has median 0.19, and none exceeds 0.8. Within
     5–10 Mb, r is 0.56–0.82 (results/split/cn_leakage_summary.json; data review #4).
   - Cost: with only 14 TEST arms, the CI is 1.33–11.2 and the primary p = 0.068.
   - Arose at: TP3, TP8.
6. **On hub-dominated benchmarks, require a gain in both folds plus a paired lower bound > 0.**
   - Evidence: focal (+0.68; folds −0.09 / +1.52) and I4 (+0.21; −0.04 / +0.51) were rejected. Focal was not
     better on TEST (Δ +0.16, CI −0.37 to 0.69).
   - The adopted G1 gain replicated on TEST.
   - Arose at: TP6, TP8.
7. **Compute cheap label-free baselines before freezing a design that must beat them.**
   - Evidence: B5 was specified at round 1, but its value (19.7, against M1's 1.59) appeared only after the seal.
     The frozen requirement "beat B5" was never met (TRAIN Δ −1.37).
   - Arose at: TP5.
8. **A within-regulator sign-permutation null treats targets as exchangeable, so it may be anti-conservative when
   predicted signs are correlated across targets.** **[not shown: an observation, not tested]**
   - Evidence: the C2 control (context fusion on permuted CN, so its signs should carry no information) scored
     concordance 0.585 against c₀ 0.514, sign p = 0.001 on HCT116 TEST (n = 354). The NC1 and C1 controls did not
     (p 0.52, 1.0).
   - REPORT.md does not discuss this.
   - Arose at: TP8.
9. **With near-collinear candidates (r ≥ 0.95), single-effect fine-mapping is decided by the prior, not the data.**
   - Evidence: in A3, σ₀ = 0.05 and df ≈ 1,335 give log ABF ≈ 0.385 z², a near-hard argmax.
   - α has within-regulator AUROC 0.498. A3's +0.36 gain came through the exp(t_cis/2) prior, i.e. regulator choice
     (M5a_ideas §3.4).
   - The algorithms reviewer predicted this at round 1 (#9).
   - Arose at: TP4, TP6.
10. **A credible single test needs physical sealing, a committed selection, a committed LOCK and a mock run on the
    same code path. Flags in code are not enough.**
    - Evidence: the sequence was seal `57f2de3`, selection `62e67b9`, LOCK `603a7f5`, run `e9085c7`.
    - The mock-plus-fresh pre-final review found 3 fixes that would have been impossible after the run.
    - Arose at: TP3, TP8.
11. **Every summary document needs an independent number-by-number check.**
    - Evidence: the report verifier checked about 120 numbers and found 4 errors and 3 over-claims, including Fisher
      p 2×10⁻³⁶ shown next to the pre-registered p 0.068.
    - The research tree had copied them.
    - Arose at: TP8.

## 3. Engineering details for reproduction

1. **Core engine: FWL residualisation per regulator chromosome.** The method is 23 pivoted-QR float64 residualisations
   plus closed-form partial correlations with a per-target cis basis (k = 2).
   - It scores 2.2×10⁸ pairs in 67 s on a 10-core laptop CPU, with a peak footprint of 6.4–8.4 GB (macOS ru_maxrss
     reports 14–20 GB).
   - It matches statsmodels OLS to 1e-8 on 480 real pairs.
   - Arose at: TP1.
2. **Guard NaNs in the ABF.** Genes with CN but no expression column had df_co = NaN. That NaN went through the log
   ABF and the segment logsumexp and wiped out α for 28.5% of all pairs. Fix: z_co = 0 for CN-only genes (D5).
   Arose at: TP6.
3. **G1 settings.** At τ = 100, 5,403 of 12,141 regulators pass. Weights w = 0 and w = 0.5 give identical top-0.1%
   results, because gated regulators fill K. Arose at: TP6.
4. **A3 attribution behaves like an argmax.** 74% of G1 trans pairs have α < 0.5, but every top-K pair has
   α ≥ 0.5. Arose at: TP6.
5. **PIE fields.**
   - `Ueffect` is a signed rank-biserial statistic in [−1, 1]. Its sign matches LFC in 99.99% of positives.
     \|Ueffect\| ≥ 0.2 keeps 47% of positives.
   - LFC magnitude is a pseudocount artefact: the perturbed mean is exactly 0 in 83.5% of expressed targets. Use
     sign and p only.
   - Arose at: TP2, TP3.
6. **Sample set and covariates.**
   - Use S2: 1,436 cell lines after removing organoids, engineered models and duplicate patients. The PureCN-complete
     S1 silently drops HCT116.
   - PureCN ploidy/WGD (24Q4) is missing for 30% of S2 (437 lines); impute and add an indicator.
   - Sequencing platform is 93% missing, so it was dropped. RNA `Stranded` was kept.
   - TP53 comes from the 24Q4 figshare matrices (Damaging > 0 or Hotspot > 0) and covers 1,255 of 1,436 lines.
   - Arose at: TP3.
7. **X chromosome.** Sex-centre X CN before applying the variable-CN rule: male X sits at log2 ≈ −1, and 95% of X
   genes pass the raw rule. Even after centring, all 716 X genes pass; the frozen flag was kept unchanged (D3).
   Arose at: TP3.
8. **Data release.** The DepMap release is *inferred* to be 26Q3. It was not md5-verified, because the portal sits
   behind a bot wall (DATA.md). Arose at: TP2.
9. **Resampling conventions.**
   - Fix the top-K threshold at its full-split value in every jackknife, bootstrap and permutation replicate.
   - Compute replicate ORs from 2×2 counts per (regulator arm × target 5 Mb bin).
   - Permute regulator labels only among regulators of the same split and line.
   - With n_perm = 200, p floors at 0.005; the final run used 1,000.
   - Arose at: TP3.
10. **Percentiles for fusion scores.** Compute them over all 18,276 DepMap-measured targets, never within the PIE
    universe: universe membership is derived from benchmark features. Arose at: TP7.
11. **Calibration.** λ_GC ≈ 2.1 on real data reflects genuine arm structure; it is 0.9–1.12 under the SIM null. The
    resulting call set (10.7 M pairs, TEST OR 1.15) is not useful, so rank instead of calling. Arose at: TP6, TP8.
12. **Reporting ORs.** Report empty-cell ORs as NA rather than Haldane-corrected values (TF 6.7 and kinase 42.9 with
    0 hits). B3, a target-agnostic score, degenerates under top-K. Arose at: TP8.
13. **Score files and evaluation.** Fused score files must carry `pair_score` and `target_weight` so that C1 can
    shuffle the pair term exactly (D6a). Evaluation uses sparse positives plus per-regulator histograms and takes
    1.4–3.2 s per score file. Arose at: TP7.
14. **Seeds.** Seeds are [20261027, 1, 2, 3, 4]. Permutation and NC1 streams come from
    `SeedSequence(20261027).spawn(n)`, with the spawn index logged. The split uses 20261027 only. Arose at: TP3.
15. **Platform and harness quirks.**
    - macOS Accelerate emits spurious matmul warnings; never promote them to errors.
    - In this agent harness, subagents cannot write report files, so the PI must save their findings (and should say
      so).
    - Arose at: TP1, TP8.

---

## Inconsistencies found between documents

1. **RESEARCH_TREE M4** still marks "[!] M1 has direction information: concordance 0.581 vs c₀ 0.476". The
   within-regulator sign permutation for that run gives p = 0.66 (null mean 0.578; `results/m4/mvp_full_hct116.txt`),
   so there is no evidence beyond each regulator's direction bias. This over-claim survived the M7 correction.
2. **RESEARCH_TREE M5a** reports "regulator-label permutation p = 0.005". That value is the 1/201 floor at
   n_perm = 200, which the ideas review flagged (§3.5). The M5 block says "≤ 0.005"; the M5a block was not reworded.
3. **DESIGN §13 (round 1)** disposes of prior-work #5 as "covered by B3/G1". The proposal also contained a
   target-expression weight, which B3 and G1 do not cover. That half returned post hoc as I1.
4. **reviews/retro/M5a_retro.md** says "~1.5 h of review rounds". The commits show the three design-review rounds
   took 22 min (18:33 → 18:55).
5. **Terminology.**
   - RESEARCH_TREE M1, DESIGN §6 and DATA_PROFILE §11 say "10% of targets give > 90% of hits" and "83% of hits come
     from non-expressed targets". There, "target" is PIE's word for the *perturbed* gene, which the project's own
     convention calls the regulator.
   - REPORT §1.6 gets this right ("perturbations").
6. **REPORT §6** contains a duplicated parenthetical: "(B4e alone scores 3.8), (B4e alone scores 3.8 on TRAIN and 4.1
   on TEST)".
7. **REPORT.md** does not mention that the C2 negative control passes the within-regulator sign test on TEST
   (p = 0.001; `results/final/summary.md`, sign table). This bears directly on the headline sign result (p = 0.012).
8. **RESEARCH_TREE M3** sub-items are still marked [~] (in progress) although M3 is [x].
9. **RESEARCH_TREE M1** cites "eQTLGen2 cis→trans vs HCT116 Perturb-seq OR ≈ 12" without the LITERATURE.md caveat:
   41 pairs, and a degree-matched null rather than a universe-based one.
