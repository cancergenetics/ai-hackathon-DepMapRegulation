# Research tree — DepMap CN → trans expression (DepMapTR)

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-27 — M2 design frozen (v3.1); M3 implementation started
The design went through three review rounds: 7 reviewers, then 5 (continued), then 1 fresh external PI. Every proposal has a disposition in DESIGN §13.
Key changes from v1:
- S2 (all 1,436 deduplicated cell lines) is the primary sample set, since S1 would drop HCT116 itself.
- The hold-out unit is the chromosome arm.
- Dependence-aware inference: regulator-label permutation, arm jackknife.
- The predicted knockdown sign uses a shrunk cis sign.
- Covariates: TP53 (24Q4), RNA batch, LOCO autosomal burden; X is sex-centred.
- Gates: G1 (cis instrument) and G2 (context expression, label-free).
- SuSiE L=1 co-segment attribution.
- Greedy tiered tuning with a strict adoption rule; TEST labels sealed outside the tree.

Running now, in parallel:
- data agent: TP53, RNA technical covariates, S2 panel + frozen universe flags (T1–T4);
- harness agent: benchmark tables, split, metrics, seal/final scripts;
- model agent: FWL engine, covariates, gates, attribution, SIM, baselines.

Interfaces are fixed in docs/INTERFACES.md. Next: review both components independently, freeze and seal, then the MVP run of M1 + baselines on TRAIN.

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
- D3 (2026-10-27): Sequencing platform dropped from Z (93% missing in S2); `Stranded` kept, though it mostly tracks sequencing era. The X variable-CN rule is kept as frozen, even though all 716 X genes pass after sex-centring (whole-X/Xi loss is common). Changing a frozen flag post hoc would be ad hoc, and G1 handles X regulators whose CN doesn't move their expression.
- D2 (2026-10-27): DESIGN v3.1 frozen. Primary = pooled top-0.1% OR on HCT116 TEST, with significance from regulator-label permutation. Tuning objective = pooled OR on TRAIN with a CMH guard. G2 tuning on pooled OR is kept (see §13, round 2).
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
- [x] M2 Design + review (3 rounds)
  - [x] DESIGN.md v1
  - [x] Review round 1 (7 reviewers → reviews/round1/) → v2
  - [x] Review round 2 (5 reviewers) → v3
  - [x] Review round 3 (fresh external PI) → v3.1 FROZEN
- [~] M3 Benchmark harness + held-out split (frozen before tuning)
  - [~] T1–T4 data tasks (TP53, RNA tech, S2 panel, frozen flags)
  - [~] benchmark.py, split, metrics, seal/final scripts
  - [~] model engine (assoc/covariates/gates/attribution/SIM), baselines
  - [ ] independent review of harness + engine
  - [ ] freeze hashes + seal TEST
- [ ] M4 Baseline models
- [ ] M5 Method development (tuning on train regulators only)
- [ ] M6 Final held-out evaluation
- [ ] M7 REPORT.md
