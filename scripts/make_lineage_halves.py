#!/usr/bin/env python
"""Write the fixed I6 lineage-half split (PI decision D6): results/split/lineage_halves.tsv.

Seeded greedy balance on n over the 32 Z lineage levels of the S2 panel (covariates.lineage_halves,
seed 20261027). Label-free; must be committed before any I6 evaluation.
"""
from pathlib import Path

import pandas as pd

from depmaptr.covariates import lineage_halves, lineage_levels
from depmaptr.io import PROC, ROOT

panel = pd.read_csv(PROC / "panel_s2.tsv", sep="\t", low_memory=False)
split = lineage_halves(panel)
n = lineage_levels(panel).value_counts()
tab = pd.DataFrame({"lineage_level": split.index, "n": n.reindex(split.index).to_numpy(), "half": split.to_numpy()})
tab = tab.sort_values(["half", "n"], ascending=[True, False])
out = ROOT / "results" / "split" / "lineage_halves.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
tab.to_csv(out, sep="\t", index=False)
print(tab.groupby("half")["n"].agg(["count", "sum"]))
print(f"wrote {out}")
