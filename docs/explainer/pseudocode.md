# DepMapTR final method — pseudocode

Selected CN-only method `g1_t100_w0__a3_r095_s2_w0` (M1 LOCO reduced form + G1 + A3), checked against
`src/depmaptr/`. n = 1,436 lines; regulators A = 12,141 variable-CN genes; targets B = 18,276 measured genes.

```text
 1  x[g] ← log2(clip(relCN[g], 1/4, 8))  for every CN gene g;   y[g] ← log2(TPM+1)
 2  Z ← lineage(32 lvls), ploidy, WGD, MSI, TP53, sex, medium, serum-free, growth, stranded (+ missing flags)
 3  arm[k] ← median x over genes on arm k  (41 arms; Xp/Xq sex-centred)
 4  for each regulator chromosome c:                                  # 23 groups (LOCO)
 5      D ← [1, Z, burden(arms ∉ c), arm[k] for k ∉ c]                # 98 columns
 6      Q ← pivoted-QR basis of D;  R(v) := v − Q Qᵀ v                # residualise once per group
 7      Ỹ ← R(y[B]) for all targets;  W̃_B ← R([x_B, 1(relCN_B<0.25) if ≥3 lines])
 8      for each CN gene g on c (all genes, needed by A3):
 9          t[g,B] ← partial t of R(x_g) vs Ỹ_B given W̃_B             # = OLS t of x_g in y_B~D+W_B+x_g
10          df[g,B] ← n − rank(D) − #active W columns − 1
11          t_cis[g] ← partial t of R(x_g) vs R(y_g)  (0 if g has no expression)
12      seg ← greedy_segments(R(x of all CN genes on c), ρ=0.95)     # join if r(gene, running mean) ≥ 0.95
13      for each segment s with m > 1 genes, each target B:
14          logpost[g] ← t_cis[g]/2 + ½log(1−λ) + ½λ·t[g,B]²,  λ = 0.05²/(0.05² + 1/df[g,B])
15          α[g,B] ← softmax_g(logpost)                               # SuSiE L=1; singleton α = 1
16  for each regulator A:
17      w_G1[A] ← 1 if t_cis[A] > 0 and t_cis[A]² ≥ 100 else 0        # G1: first-stage F ≥ 100
18      s_A ← −1 if t_cis[A] ≤ −3 else +1
19      for each target B on another chromosome:
20          S[A,B] ← t[A,B] · w_G1[A] · α[A,B]                         # final score; rank by |S|
21          kd_sign[A,B] ← −sign(S[A,B]) · s_A                         # predicted change of B on knockdown of A
22  # call set (reported, not used for ranking): one lineage-stratified permutation of R(x) rows (NC1);
23  t* ← min threshold with  #(|t_NC1| ≥ t*) / #(|t| ≥ t*) ≤ 0.05;   call[A,B] ← |t[A,B]| ≥ t*
24  # optional context layer (reported variant I1, γ=2; not the selected method):
25  ρ_P[A,B] ← percentile of |S[A,B]| among all finite trans entries (0 where S = 0)
26  ρ_T[B]   ← percentile of B's DepMap HCT116 log2(TPM+1) among all targets
27  S_ctx[A,B] ← sign(t[A,B]) · ρ_P[A,B] · ρ_T[B]²
28  # evaluation (benchmark.py): top 0.1% of trans universe pairs by |S| vs Perturb-seq p_adj < 0.05
```

## Step → code map (all names verified to exist)

| pseudocode line(s) | step | file:function |
|---|---|---|
| 1 | CN transform | `io.py:cn_transform`; data loaded by `io.py:load_real` |
| 2 | shared covariates Z | `covariates.py:build_Z`, `covariates.py:lineage_levels` |
| 3 | arm-level CN | `covariates.py:arm_cn` (`sex_centre` for X) |
| 5 | LOCO design D_c | `covariates.py:DesignBuilder.design`, `covariates.py:burden` |
| 6 | residualisation | `assoc.py:Residualiser` (pivoted QR, rank/cond checks) |
| 7 | residualised targets and cis basis | `pipeline.py:prepare` (XT, YT, DEL), `pipeline.py:GroupCache.get` |
| 9–10 | reduced-form t and df (FWL sweep) | `pipeline.py:group_t` → `assoc.py:chunked_pair_partial` → `assoc.py:pair_partial` |
| 11 | cis t for G1 | `assoc.py:diag_partial` (called in `pipeline.py:run`; also `pipeline.py:cis_t`) |
| 12 | co-segments | `attribution.py:greedy_segments` |
| 13–15 | SuSiE L=1 attribution α | `attribution.py:segment_alpha`, `attribution.py:log_abf` |
| 17 | G1 gate | `gates.py:g1_pass`, `gates.py:g1_weight` |
| 18 | knockdown sign rule s_A | `gates.py:s_cis` |
| 20 | final score | `pipeline.py:compose_score` (orchestrated by `pipeline.py:run`) |
| 21 | predicted knockdown sign | `benchmark.py:Evaluation.pred_kd_sign` (`gates.py:predicted_kd_sign` for tests) |
| 22 | NC1 permutation | `pipeline.py:nc1_permutation` → `assoc.py:lineage_permutation` |
| 23 | call set | `calls.py:trans_mask`, `calls.py:t_star`, `calls.py:call_matrix` (driven by `scripts/run_model.py:main`) |
| 25–27 | context layer | `pipeline.py:context_fusion`, `pipeline.py:trans_percentile`, `pipeline.py:target_context_percentile` |
| 28 | evaluation | `benchmark.py:Evaluation.topk`, `benchmark.py:compute_metrics`, `benchmark.py:evaluate` |
| — | CLI / variant composition | `scripts/run_model.py:resolve`, `results/variants.yaml` (`g1_t100_w0` + `a3_r095_s2_w0`) |
