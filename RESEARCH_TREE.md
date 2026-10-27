# Research tree — DepMap CN → trans expression (DepMapTR)

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-27 — M0 project setup
Project initialised: PROMPT.md saved verbatim, git repo, Python venv (.venv, packages locked in
requirements.lock.txt). Compute: 10 CPU, 32 GB RAM, no GPU. Starting data: DepMap omics (WGS gene CN,
protein-coding expression log2(TPM+1)), Model/ModelCondition/OmicsProfiles/SubtypeTree metadata, HGNC protein-coding list.
Benchmark: PIE X-Atlas/Orion DE tables (HCT116, HEK293T) downloading, pinned to HF commit abfa3da.

Next: parallel literature review (LITERATURE.md) and data profiling (DATA.md, release identification,
covariate availability), then DESIGN.md and up to three review rounds.

## Decision notes
- D0 (2026-10-27): Use pinned HF commit abfa3da0f3bb98c9df9c90082386ed0c7c9282d3 for the benchmark; use only `de/*.parquet`
  (pre-computed perturbation–mRNA associations), never raw Perturb-seq counts.

## Tree
- [~] M0 Setup
  - [x] PROMPT.md, git, venv
  - [~] Download benchmark DE tables
- [ ] M1 Literature + data profiling
  - [ ] LITERATURE.md
  - [ ] DATA.md (source, version, checksum, license)
  - [ ] Data profile: sample overlap, covariates, benchmark structure
- [ ] M2 Design + review (≤3 rounds)
- [ ] M3 Benchmark harness + held-out split (frozen before tuning)
- [ ] M4 Baseline models
- [ ] M5 Method development (tuning on train regulators only)
- [ ] M6 Final held-out evaluation
- [ ] M7 REPORT.md
