# Research tree — DepMap CN → trans expression (DepMapTR)

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-27 — M5a: pre-registered tiers complete (23 of 30 configurations used)
Greedy tuning on HCT116 TRAIN adopted G1 (cis instrument gate, F ≥ 100) and A3 (SuSiE L=1 co-segment attribution, ρ=0.95, s=2, no co-expression). The incumbent `g1_t100_w0__a3_r095_s2_w0`:
- HCT116 TRAIN: OR@0.1% = 5.0 (M1 1.59), regulator-label permutation p = 0.005, degree-stratified permutation p = 0.055, CMH 1.09.
- Sign concordance 0.561 vs c0 0.483 (within-regulator permutation p = 0.005); sign-concordant OR 6.2.
- HEK293T TRAIN: OR 2.5.

Rejected: arm-handling alternatives (unstable), G2 (redundant after G1), target standardisation, hinge, CN floor, RINT.

The method is still far below the DepMap-only hub baseline B5 (regulator cis×CN-SD × target HCT116 expression: OR 19.7). It also falls short on CMH (Δlog −0.93). Most pooled enrichment in this benchmark comes from regulator hubness and target detectability. The CN pair evidence adds little within-regulator target discrimination (CMH ≈ 1).
Next (D6): a small post-hoc tier that fuses pair evidence with DepMap-only target-context priors, testing whether CN pair evidence adds anything beyond B5. Also extend G1 τ past the grid maximum. A fresh agent runs the retrospective and a fresh scientist proposes ideas.

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
- D5 (2026-10-27, before any evaluation): code review found two problems in A3 attribution. (a) CN-only genes produced NaN co-expression evidence that spread across whole segments; they now get z_co = 0. (b) G2 was lost for singletons under A3. The final score is now always S = t·w_G1·w_G2·α, with π_g ∝ exp(t_cis/s)·w_G2(g) inside A3 (a design fix of §5.6). The SIM is extended so that each confounder test fails when its adjustment is removed. All pre-fix score runs are deleted.
- D4 (2026-10-27, before any model result): PIE `Ueffect` is a signed rank-biserial effect in [−1, 1] whose sign matches LFC in 99.99% of positives. The v3.1 sensitivity rule |Ueffect − 0.5| ≥ 0.1 assumed a [0, 1] U scale and kept 98.6% of positives, so it was meaningless. It was translated to the equivalent |Ueffect| ≥ 0.2, which keeps 47% of HCT116 TRAIN positives. Sensitivity set only; the primary label is unchanged.
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
- [x] M3 Benchmark harness + held-out split (frozen before tuning)
  - [~] T1–T4 data tasks (TP53, RNA tech, S2 panel, frozen flags)
  - [~] benchmark.py, split, metrics, seal/final scripts
  - [~] model engine (assoc/covariates/gates/attribution/SIM), baselines
  - [x] independent review of harness (no blockers) + engine (1 blocker fixed, D5)
  - [x] freeze hashes + seal TEST (2026-10-27; ~/DepMapTR_sealed, results/FROZEN.sha256)
- [x] M4 Baseline models (HCT116 TRAIN, results/m4/)
  - [!] M1 pooled OR@0.1% = 1.59 (permutation p 0.25; CMH 1.10). Hub baselines dominate: B5 (reg cis×CN SD × target HCT116 expression) 19.7, B4e (target HCT116 expression) 7.1, B2 (co-expression) 5.8. Controls OK: B0 0.98, NC1 0.75, B6 (arm-only) 0.97.
  - [!] M1 has direction information: concordance 0.581 vs independence c0 0.476.
  - [x] B4e CMH = 3.8: the expression-decile strata do not fully remove target detectability, so CMH is not a pure pair-specific measure.
- [~] M5 Method development (tuning on train regulators only; greedy tiers, DESIGN §5.5, results/tuning/)
  - [-] Tier 1 arm handling (focal, loao, cnpcs20): none adopted. Focal Δlog-OR +0.68, but fold1 −0.09 / fold2 +1.52, jk lb90 −0.62, so the gain is unstable. LOCO retained.
  - [!] Tier 2 G1 instrument gate: all adopted under §8.2. τ=100 is best: Δlog-OR +0.79 (folds 0.65/1.12, jk lb90 +0.61), ΔCMH −0.05. Adopted g1_t100_w0 (w=0.5 identical to 2 d.p.). The gain is monotone in τ up to the grid maximum.
  - [-] Tier 3 G2 context gate on top of G1: not adopted (Δlog-OR +0.05, lb90 ≈ 0). Regulators that pass a strong cis gate are already expressed. Incumbent g1_t100_w0: OR 3.50 (jk CI 2.5–4.8), perm p 0.03, CMH 1.04, concordance 0.578 vs c0 0.487.
  - [!] The G1 gain is regulator selection, not within-regulator target ranking (CMH ≈ 1).
  - [x] Tier 4 A3 attribution: a3_r095_s2_w0 adopted (Δlog-OR +0.36, folds 0.43/0.28, lb90 +0.04). w_co=0.25 also passes with a smaller Δ; ρ=0.98 and s=5 do not.
  - [-] Tier 5 target standardisation: Δ −0.03, not adopted
  - [-] Tier 6 hinge (loss slope): Δ −0.11, not adopted
  - [-] Tier 7 floor 2^-5 (Δ −0.01), RINT (Δ −0.16): not adopted
  - [!] Incumbent after pre-registered tiers: g1_t100_w0__a3_r095_s2_w0, HCT116 TRAIN OR 5.0 (perm p 0.005), sign concordance significant; still well below B5 (19.7)
- [ ] M6 Final held-out evaluation
- [ ] M7 REPORT.md
