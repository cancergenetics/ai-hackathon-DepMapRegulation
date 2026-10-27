# Research tree — DepMap CN → trans expression (DepMapTR)

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-27 — M1 literature + data done; M2 design review round 1 in progress
LITERATURE.md (59 verified refs) shows that observational trans sets enrich weakly in Perturb-seq (OR ~1–2), except
small causally anchored sets (eQTLGen2 cis→trans vs HCT116 Perturb-seq OR ≈ 12). The main levers: don't regress out expression PCs,
use LOCO arm CN for cross-chromosome confounding, gate regulators on the cis CN→mRNA instrument, attribute
co-segmented genes, and control hub/degree structure in the evaluation. The data foundation is built and tested (8 tests); the release is inferred as DepMap 26Q3.
There are 1,724 matched models (1,509 cell lines); PureCN ploidy/WGD (24Q4) cover 1,037. CN is linear relative CN, ploidy-normalised.

The benchmark is hub-dominated: 10% of targets give >90% of trans hits. In HEK293T, 83% of hits come from targets that are not expressed
(likely artefacts), so HCT116 is primary. DESIGN.md v1 was drafted with a frozen benchmark definition, a 5 Mb-bin regulator hold-out, and the FWL partial-correlation
core model with LOCO arms. Seven independent reviewers (external PI, statistics, data audit, domain biology, prior work,
algorithms/compute, AI failure modes) are reviewing it now.

### 2026-10-27 — M0 project setup
Project initialised: PROMPT.md saved verbatim, git repo, Python venv (.venv, packages locked in
requirements.lock.txt). Compute: 10 CPU, 32 GB RAM, no GPU. Starting data: DepMap omics (WGS gene CN,
protein-coding expression log2(TPM+1)), Model/ModelCondition/OmicsProfiles/SubtypeTree metadata, HGNC protein-coding list.
Benchmark: PIE X-Atlas/Orion DE tables (HCT116, HEK293T) downloading, pinned to HF commit abfa3da.

Next: parallel literature review (LITERATURE.md) and data profiling (DATA.md, release identification,
covariate availability), then DESIGN.md and up to three review rounds.

## Decision notes
- D1 (2026-10-27): HCT116 is the primary benchmark line; HEK293T is secondary (non-DepMap line, artefact-dominated hits). Decided from data profiling, before any model results.
- D0 (2026-10-27): Use pinned HF commit abfa3da0f3bb98c9df9c90082386ed0c7c9282d3 for the benchmark; use only `de/*.parquet`
  (pre-computed perturbation–mRNA associations), never raw Perturb-seq counts.

## Tree
- [x] M0 Setup
  - [x] PROMPT.md, git, venv
  - [x] Download benchmark DE tables (sizes verified)
- [x] M1 Literature + data profiling
  - [x] LITERATURE.md
  - [x] DATA.md (source, version, checksum, license); release inferred 26Q3 (not md5-verified: portal bot-wall)
  - [x] Data profile (docs/DATA_PROFILE.md)
  - [!] HEK293T benchmark: 83% of trans hits come from non-expressed targets → HCT116 primary; non-expressed targets used as a null
  - [!] Benchmark hub-dominated: top 10% of targets give >90% of hits → report CMH-stratified OR + hub baselines
- [~] M2 Design + review (≤3 rounds)
  - [x] DESIGN.md v1
  - [~] Review round 1 (7 reviewers → reviews/round1/)
- [ ] M3 Benchmark harness + held-out split (frozen before tuning)
- [ ] M4 Baseline models
- [ ] M5 Method development (tuning on train regulators only)
- [ ] M6 Final held-out evaluation
- [ ] M7 REPORT.md
