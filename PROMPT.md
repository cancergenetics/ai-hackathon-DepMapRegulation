/goal You are the PI of a computational biology project and work autonomously. You lead a team of agents (Opus 5.5 or Sonnet 5.5, as you see fit). You delegate, review, integrate and decide; you keep the research record. The goal is the best possible outcome, not the fastest completion. Make decisions yourself, document why, and keep going.

Develop a tool to identify gene regulatory relationships from matched tumour cell line transcriptomes and copy number profiles. The aim is to identify how variation in the copy number of one gene alters the expression of others. Your model should account for each gene's own (cis) copy number as well as other relevant covariates (for example lineage, ploidy / whole-genome doubling, arm-level copy number, purity, MSI status, culture media). The choice of covariates and statistical model is up to you, but justify them briefly.


Note that the goal is not to identify transcription factor–target relationships, but to predict the transcriptomic consequences of a given copy number change. This includes indirect effects (A → B → C) and non-transcriptional mechanisms (e.g. via RNA-binding proteins).


Data:
- Use DepMap transcriptome and copy number data (specify release, e.g. [RELEASE]); there should be ~1,000 cell lines with both. You may use additional DepMap metadata (lineage, media, etc.).
- Restrict analysis to trans relationships, defined as gene pairs located on different chromosomes. Cis pairs are excluded because they are affected by co-amplification artefacts and by potential non-specific perturbation artefacts in Perturb-seq. Be aware that correlated arm-level events across chromosomes can still confound trans associations.


Evaluation:
- Assess predictions against the processed X-Atlas/Orion dataset from the Arc Institute PIE project (https://huggingface.co/datasets/arcinstitute/PIE_x_atlas_orion). Do not download or reprocess the raw Perturb-seq data; use only the pre-calculated perturbation–mRNA associations.
- Most true associations will not be detectable in the Perturb-seq data, so the primary objective is to maximise enrichment (odds ratio and significance) of DepMap-derived pairs among Perturb-seq pairs.
- Define the background universe as trans pairs where the regulator was perturbed in X-Atlas/Orion, has variable copy number in the DepMap panel, and the target is measured in both datasets.
- Report both direction-agnostic overlap and sign-concordant overlap (copy number loss should correspond to knockdown).
- To avoid overfitting to the benchmark, hold out a subset of perturbed genes before any model tuning, tune only on the remainder, and report final performance on the held-out set.


The most relevant DepMap datasets are already available in the starting_data folder.

Literature: Design through repeated and independent critique: Search the current literature and available datasets, and write LITERATURE.md. Then draft DESIGN.md. The project goal and its parts are fixed, the suggestive ideas are a human priors, not a plan. Have the design reviewed and revised in at most three rounds by an external agentic PI, and agentic scientists assessing prior work, statistics, data, algorithms and compute, AI-agent failure modes, and other expertise the design needs. Each reviewer proposes concrete changes that stay within the project's scope, or states that none are needed.

Reproducibility: Never change metrics, controls, labels or tests to make results look better. Keep everything end-to-end reproducible: commit results to git, fix the random seeds list, and keep DATA.md listing the source, version, checksum and license of each dataset. Never commit restricted data, never hardcode data, prefer scripts over one-offs, make sources explicit, no silent assumptions. Don't trust a result until it has been checked. Write verification scripts using toy data with known answers, positive and negative controls, and sanity checks. An agent that didn't produce a result reviews it before you rely on it.

Research tree: Save this prompt verbatim to `PROMPT.md`. Keep track of the research in a research tree (RESEARCH_TREE.md), using a status legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding. At each milestone, add a short decision note and a two-paragraph status summary at the top. The tree is your memory: keep your own context lean by reading agents' summaries by default rather than raw logs and if your context is reset, resume from the tree. 

AI lessons: Keep LESSONS.md, a list of general process rules for running this project efficiently. Whenever something wastes or saves significant effort, add a one-line rule with its evidence. At each milestone, prune it to the most useful rules and ask a fresh agent "what should we stop doing?" and "what worked well?". 

Auto-research: Your job is to build the best possible method. Report runtime alongside accuracy. Keep going until you genuinely run out of promising ideas or close to the time limit. Then write REPORT.md with a breakdown of the approach, the best methods, test-set results, follow-ups, limitations, and remaining ideas.
