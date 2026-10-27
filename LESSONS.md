# LESSONS — process rules for running this project efficiently

Each rule is one line with its evidence. Pruned at each milestone.

1. Check the environment (CPU, RAM, packages, network) before planning compute. Evidence: system Python had no numpy; a venv was needed first (M0).
2. Read the HF dataset README and file sizes before downloading. Evidence: the DE tables are 5.2 GB, and the preprocessed arrays (unneeded) would have added another 6 GB (M0).
