# What was learned

> AI-generated project, not yet verified by a human. See the Disclaimer in the [README](../../README.md).
> Method: [METHOD.md](METHOD.md) · Journey: [JOURNEY.md](JOURNEY.md) · Interactive version: [index.html](index.html)

**The main lesson: copy-number associations, if they say anything here, say which regulators matter rather than which
targets respond (weakly supported), and a pooled enrichment score on a single-line Perturb-seq screen mostly rewards
hubs and detectable targets.** Insights are ranked
within each category by how widely they apply, then by surprise. **[not shown]** marks a belief this project did not
demonstrate. Each insight links to the [turning point](JOURNEY.md#turning-points) where it arose.

## The top two in each category

**Scientific**
1. **Which regulators, not which targets.** Regulators whose own mRNA strongly tracks their copy number are more
   likely to have detectable trans effects on knockdown (gate: +0.97 log-OR on held-out arms, CI 0.30–1.64; the overall
   enrichment itself was not significant). Among a regulator's targets, copy number ranks responders at chance
   (AUROC 0.500). [More](#sci-1)
2. **A regulator's trans row is mostly its arm's signature.** A gated regulator's t-statistics correlate with its arm's
   at median r = 0.87, so neighbours on an arm share nearly the same predictions. [More](#sci-3)

**Methodological**
1. **Split pooled enrichment into its parts before crediting a method:** regulator choice, target detectability and
   pair-specific evidence. A +1.57 log-OR gain here contained no detectable pair-level knowledge. [More](#meth-1)
2. **Tune on a metric that measures the claim.** Tuning on pooled odds ratio adopted only regulator-level gains and
   rejected the one option that improved target ranking. [More](#meth-2)

## 1 · Scientific and domain insights

1. <a id="sci-1"></a>**DepMap copy number carries no detectable information about which targets of a given regulator respond to its
   knockdown.** AUROC 0.507 (training, 333 regulators) and 0.500 (test, 179); stratified odds ratio 0.94 on test (CI
   0.77–1.15, any lift below ~15%); target shuffle 3.76 vs 3.85; adjusting A's own arm does no better (0.511).
   *Arose at [TP7](JOURNEY.md#tp7), confirmed at [TP8](JOURNEY.md#tp8).*
2. <a id="sci-2"></a>**A regulator's cis dosage response predicts whether its knockdown has detectable trans consequences, even though
   it is measured across all lineages.** Gate (F ≥ 100) +0.79 log-OR on training, rising with the threshold; +0.97 (CI
   0.30–1.64) on test; target-level effect −0.05. *Arose at [TP4](JOURNEY.md#tp4), [TP6](JOURNEY.md#tp6); replicated
   at [TP8](JOURNEY.md#tp8).*
3. <a id="sci-3"></a>**A regulator's trans-association row is mostly the signature of its chromosome arm.** Median r = 0.87 with the
   arm-only row; 91% of rows exceed 0.7; 54 regulators on 1q share 1,308 top-ranked targets. Arm copy number alone does
   not enrich (training 0.97, test 0.61). *Arose at [TP7](JOURNEY.md#tp7).*
4. <a id="sci-4"></a>**In a single-line CRISPRi Perturb-seq screen, regulator hubness and target expression depth largely decide which
   pairs count as hits.** A label-free hub prior scores 19.7 / 22.8 (training / test) against 5.00 / 3.85 for the
   copy-number method; even the fusion built on shuffled copy number (10.6) beats it on test (−1.01, CI −1.99 to
   −0.04). *Arose at [TP2](JOURNEY.md#tp2), [TP5](JOURNEY.md#tp5); quantified at [TP7](JOURNEY.md#tp7).*
5. <a id="sci-5"></a>**Whether a regulator is expressed in the screened line changes its hit yield about 20-fold; a population-level
   instrument is not a context-level one.** Per-pair hit rate 0.0006 → 0.0124 across HCT116 expression bins of gated
   regulators; MYCN passes the gate via neuroblastoma amplicons, is not expressed in HCT116, and took 1,701 top-ranked
   training slots for 1 hit. *Arose at [TP2](JOURNEY.md#tp2), [TP7](JOURNEY.md#tp7).*
6. <a id="sci-6"></a>**Whether copy number predicts the direction of a knockdown effect is not established.** Test concordance 71.8% vs
   64.7% (p = 0.012, n = 131, 6 regulators with ≥ 5 overlaps); the shuffled-copy-number fusion control passes the same
   test (58.5% vs 53.2%, p = 0.001); n.s. in the 999-line set (p 0.13) and least-leaky arms (p 0.48). *Arose at
   [TP8](JOURNEY.md#tp8); corrected 28 Oct (decision note D9).*
7. <a id="sci-7"></a>**Trans t-statistics skew negative, more so for highly expressed targets; part is a TPM compositional artefact.**
   t ≤ −3 is twice as common as t ≥ 3 (3.5% vs 1.7%); mean t −0.075 → −0.311 across target-expression quintiles;
   renormalising TPM without A's chromosome moves the chr1 row mean from −0.12 to +0.09 but does not change the ranking.
   The remainder is attributed to an own-chromosome aneuploidy response **[not shown]**. *Arose at
   [TP7](JOURNEY.md#tp7).*
8. <a id="sci-8"></a>**HEK293T trans hits in X-Atlas/Orion are dominated by perturbations of non-expressed genes.** 83.0% vs 12.7% in
   HCT116, with a different up/down balance (59% up); the method fails there (OR 0.49, 6 overlaps). That these hits are
   artefacts is inferred **[not shown]**. *Arose at [TP2](JOURNEY.md#tp2).*
9. <a id="sci-9"></a>**The copy-number signal is not centred on transcription factors** **[not shown: low power]**. RNA-binding-protein
   regulators OR 3.56 (30 top-ranked hits); TFs 0 hits (137 regulators); ribosomal 0 hits in 631 top-ranked pairs.
   *Arose at [TP8](JOURNEY.md#tp8).*

## 2 · Methodological insights that generalise

1. <a id="meth-1"></a>**Split any pooled enrichment into regulator allocation, target detectability and pair-specific evidence before
   crediting a method.** The detectability fusion gained +1.57 and passed the strict rule, but the exact
   within-regulator target shuffle gave the same OR (24.5 vs 24.1); ~94% of the best model's log-enrichment was
   regulator allocation. *Arose at [TP5](JOURNEY.md#tp5), [TP7](JOURNEY.md#tp7).*
2. <a id="meth-2"></a>**Choose the tuning objective from the scientific claim, not the headline metric.** Pooled-OR tuning adopted only
   allocation gains (gate: target-level −0.05) and rejected the only option that improved target ranking (+0.42
   target-level, pooled lower bound below 0). A reviewer had proposed a target-level objective at round 1. *Arose at
   [TP3](JOURNEY.md#tp3), [TP6](JOURNEY.md#tp6), [TP7](JOURNEY.md#tp7).*
3. <a id="meth-3"></a>**Give every control, metric and threshold a known-answer test before first real use.** A loader bug made the
   target-shuffle control report a spurious pair-specific gain (+0.063; exact −0.02); a frozen effect-size filter
   assumed the wrong scale and kept 98.6% of positives; one adoption rule could not be met by construction. *Arose at
   [TP7](JOURNEY.md#tp7) (and [TP3](JOURNEY.md#tp3)).*
4. <a id="meth-4"></a>**A within-regulator sign-shuffle test treats targets as exchangeable, so it can be anti-conservative when targets
   respond in correlated modules.** The shuffled-copy-number fusion control scored p = 0.001 (n = 354); the plain
   shuffled-copy-number control did not (p 0.52). The module explanation is an inference **[not shown]**; a calibrated
   test would permute regulators rather than target signs. *Arose at [TP8](JOURNEY.md#tp8); decision note D9.*
5. <a id="meth-5"></a>**Stratified "pair-specific" statistics must be stratified on the context that drives detectability.** Stratifying
   by pan-DepMap expression deciles still let target expression alone score 3.8 (training) and 4.1 (test); HCT116
   deciles reduced but did not remove this (3.4). *Arose at [TP5](JOURNEY.md#tp5), [TP7](JOURNEY.md#tp7).*
6. <a id="meth-6"></a>**Hold out copy-number-derived predictions by chromosome arm, and budget for the precision this costs.** Median
   highest test–training CN correlation 0.19 (none > 0.8) vs 0.56–0.82 within 5–10 Mb; cost: 14 test arms, CI
   1.33–11.2, p = 0.068. *Arose at [TP3](JOURNEY.md#tp3), [TP8](JOURNEY.md#tp8).*
7. <a id="meth-7"></a>**On hub-dominated benchmarks, require a gain in both folds plus a paired lower bound above zero.** It rejected
   focal (+0.68, folds −0.09 / +1.52; not better on test) and coarser blocks; the gain it adopted replicated. *Arose at
   [TP6](JOURNEY.md#tp6), [TP8](JOURNEY.md#tp8).*
8. <a id="meth-8"></a>**Compute cheap label-free baselines before freezing a design that must beat them.** The hub prior was specified at
   round 1 but first computed after the seal (19.7 vs 1.59); the requirement to beat it was never met (−1.37). *Arose at
   [TP5](JOURNEY.md#tp5).*
9. <a id="meth-9"></a>**With near-collinear candidates (r ≥ 0.95), single-effect fine-mapping is decided by the prior.** With σ₀ = 0.05
   and ~1,335 df, log Bayes factor ≈ 0.385 z², a near-hard choice; attribution's within-regulator AUROC 0.498. *Arose at
   [TP4](JOURNEY.md#tp4), [TP6](JOURNEY.md#tp6).*
10. <a id="meth-10"></a>**A credible single test needs physical sealing, a committed selection, a committed lock and a mock run on the same
    code path.** Seal `57f2de3` → selection `62e67b9` → lock `603a7f5` → run `e9085c7`; the mock and a fresh pre-final
    review found 3 fixes. *Arose at [TP3](JOURNEY.md#tp3), [TP8](JOURNEY.md#tp8).*
11. <a id="meth-11"></a>**Every summary document needs an independent number-by-number check.** ~120 numbers checked, 4 errors and 3
    over-claims found (including Fisher p 2×10⁻³⁶ next to the pre-registered 0.068); the direction over-claim survived
    until a later control check. *Arose at [TP8](JOURNEY.md#tp8).*

## 3 · Engineering details for reproduction

1. **Core engine:** one pivoted-QR residualisation per regulator chromosome (23), then closed-form partial
   correlations; 2.2×10⁸ pairs in 67 s on a 10-core laptop; matches OLS to 1e-8. *[TP1](JOURNEY.md#tp1)*
2. **Guard missing values in the Bayes-factor step:** genes with copy number but no expression produced a missing df
   that wiped out α for 28.5% of pairs; fixed by setting their co-expression z to 0. *[TP6](JOURNEY.md#tp6)*
3. **Gate weight:** weights 0 and 0.5 for failing regulators give identical top-0.1% results; 5,403 of 12,141 pass.
   *[TP6](JOURNEY.md#tp6)*
4. **PIE tables:** use sign and p only (fold changes are pseudocount artefacts); the effect statistic is a signed
   rank-biserial in [−1, 1], |value| ≥ 0.2 keeps 47% of positives. *[TP2](JOURNEY.md#tp2), [TP3](JOURNEY.md#tp3)*
5. **Samples and covariates:** use the 1,436-line set; impute missing ploidy/doubling (30%) with an indicator; drop
   sequencing platform (93% missing); TP53 covers 1,255 of 1,436 lines; sex-centre X copy number before the
   variable-copy-number rule. *[TP3](JOURNEY.md#tp3)*
6. **Resampling:** fix the top-K threshold at its full-split value in every replicate; permute regulator labels within
   split and line; 200 permutations floor p at 0.005 (final run used 1,000); seeds [20261027, 1, 2, 3, 4].
   *[TP3](JOURNEY.md#tp3)*
7. **Fusion percentiles** over all DepMap-measured targets, never within the benchmark universe; store the pair term and
   target weight separately so the target shuffle is exact. *[TP7](JOURNEY.md#tp7)*
8. **Rank rather than call:** λ ≈ 2.1 on real data reflects arm structure (0.9–1.12 under a simulated null); the 5%-FDR
   set (10.7 M pairs) has test OR 1.15. *[TP6](JOURNEY.md#tp6), [TP8](JOURNEY.md#tp8)*
9. **Report empty-cell odds ratios as not available:** Haldane-corrected 6.7 and 42.9 for classes with 0 hits are
   artefacts. *[TP8](JOURNEY.md#tp8)*
10. **Platform quirks:** macOS Accelerate emits spurious matmul warnings; in this harness sub-agents could not write
    report files, so the main agent saved condensed versions; the DepMap release could not be checksum-verified.
    *[TP1](JOURNEY.md#tp1), [TP2](JOURNEY.md#tp2), [TP8](JOURNEY.md#tp8)*
