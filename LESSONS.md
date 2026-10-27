# LESSONS — process rules for running this project efficiently

Each rule is one line with its evidence. Pruned at each milestone.

1. Check the environment (CPU, RAM, packages, network) before planning compute. Evidence: system Python had no numpy; a venv was needed first (M0).
2. Read the HF dataset README and file sizes before downloading. Evidence: the DE tables are 5.2 GB, and the preprocessed arrays (unneeded) would have added another 6 GB (M0).
3. Run independent reviewers in parallel, each with a narrow remit; issues several of them raise are almost certainly real. Evidence: 4 of 7 round-1 reviewers independently flagged pair dependence in Fisher p, and 4 flagged the 5 Mb split leaking via arm-level CN.
4. Have a second agent audit the data agent's numbers before the design relies on them. Evidence: the round-1 audit found that S1 (PureCN-complete) silently excludes HCT116, the benchmark's own cell line, and 45% of MSI-H lines (M2).
5. Continue round-1 reviewers via SendMessage for round 2 instead of spawning new ones: they check their own items cheaply. Use one fresh reviewer for the final round, to catch reviewers anchoring on their own earlier proposals.
6. Sex chromosomes break naive CN-variability and burden statistics (male X ≈ log2 −1); always sex-centre X before thresholding. Evidence: 3 round-2 reviewers caught this; 95% of X genes passed the raw rule.
