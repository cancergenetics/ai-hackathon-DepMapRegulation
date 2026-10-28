# DepMapTR — trans-transcriptomic consequences of copy-number change (AI hackathon project)

## Disclaimer

This is an exercise in assessing the bioinformatics capabilities of agentic AI. This particular exercise was run by the Ryan lab at University College Dublin as part of a larger hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects. All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication and bioinformatics tool.

Models and machine: the project was run in Claude Code. The PI agent and most sub-agents (data, benchmark harness, model engine, literature, design and code reviewers) used Claude Opus 5.5 (`claude-opus-5-5`); some reviewers and the retrospectives used Claude Sonnet 5.5 (`claude-sonnet-5-5`). All computation ran locally on an Apple M4 laptop (10 CPU cores, 32 GB RAM, no GPU).

Currently none of the results have been verified by a human. We plan to verify the agent’s methodology and benchmarking, and we will update this section accordingly.

The rest of this repository is entirely AI generated. It includes code for implementation of the predictive model and a AI-written report in the style of a journal article, detailing the final method, the path followed to reach it, and benchmarking against other methods.

## What is here

| file / directory | contents |
|---|---|
| `PROMPT.md` | the project brief, verbatim |
| `REPORT.md` | final report: approach, results on held-out regulators, limitations, follow-ups |
| `RESEARCH_TREE.md` | research log: status summaries, decision notes D0–D8a, tree of what was tried |
| `DESIGN.md` | design (v3.1, frozen after three review rounds) and dispositions of every review proposal |
| `LITERATURE.md` | literature review |
| `DATA.md` | data sources, versions, checksums and licences |
| `LESSONS.md` | process lessons for running agentic research projects |
| `docs/` | data profile, benchmark and model implementation notes, interfaces |
| `reviews/` | all design, code, pre-final and report reviews, idea and retrospective reviews |
| `src/depmaptr/` | the tool (association engine, covariates, gates, attribution, simulation, benchmark harness) |
| `scripts/` | data preparation, model runs, evaluation, final evaluation |
| `tests/` | unit, simulation (known-answer) and real-data invariant tests |
| `results/` | committed metrics: tuning log, TRAIN comparisons, split definitions, the single locked TEST evaluation (`results/final/`) |

## Data

No data files are included in this repository. To reproduce, obtain the following and follow DATA.md and REPORT.md §8:
- the DepMap files listed in DATA.md (placed in `starting_data/`);
- the Arc Institute PIE X-Atlas/Orion DE tables (Hugging Face `arcinstitute/PIE_x_atlas_orion`, commit `abfa3da`; CC BY-NC-SA 4.0).

All data and derived score matrices are git-ignored.
