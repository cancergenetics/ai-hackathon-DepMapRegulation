# Final TEST evaluation (MOCK: TRAIN fold-2 as pseudo-TEST)

## pooled

| role | run | label | line | split | subset | n_pairs | K | OR@0.1% | OR@0.01% | OR@1% | fisher_p(indep.) | jk_CI95_lo | jk_CI95_hi | perm_p | perm_p_degree | OR/median_null | OR/median_null_degree |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 4.11 | 5.73 | 5 | 6.15e-47 | 2.44 | 6.92 | 0.034 | 0.0769 | 2.27 | 1.59 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 1.7 | 2.07 | 1.87 | 4.26e-05 | 0.755 | 3.85 | 0.0939 | 0.105 | 1.73 | 1.45 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 26.2 | 32.9 | 18.3 | 0 | 19.2 | 35.7 | 0.000999 | 0.000999 | 2.2 | 1.65 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 6.99 | 6.26 | 6.5 | 6.58e-127 | 2.43 | 20.1 | 0.000999 | 0.015 | 4.85 | 3.04 |
| baseline | B0 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 1.06 | 1.03 | 0.992 | 0.383 | 0.794 | 1.41 | 0.378 | 0.505 | 1.06 | 1.01 |
| baseline | B1 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 1.26 | 2.85 | 1.44 | 0.0635 | 0.677 | 2.36 | 0.845 | 0.815 | 0.638 | 0.707 |
| baseline | B2 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 5.38 | 7.59 | 4.63 | 1.97e-79 | 3.57 | 8.12 | 0.207 | 0.533 | 1.73 | 0.956 |
| baseline | B3 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 4.64 | 0.128 | 4.17 | 8.53e-60 | 0.000525 | 40995 | 0.143 | 0.397 | 56.9 | 2.07 |
| baseline | B4 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 5.41 | 10.8 | 3.74 | 3.64e-80 | 4.58 | 6.39 | 0.474 | 0.462 | 1 | 1.01 |
| baseline | B4e |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 8.36 | 9.19 | 7.53 | 4.63e-172 | 6.87 | 10.2 | 0.441 | 0.482 | 1.01 | 1 |
| baseline | B5 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 16.4 | 28.5 | 16.1 | 0 | 6.1 | 44.3 | 0.201 | 0.696 | 1.54 | 0.832 |
| baseline | B6 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 0.902 | 0.257 | 1.07 | 0.751 | 0.381 | 2.14 | 0.346 | 0.407 | 1.07 | 1.03 |
| control | NC1 | NC1: residual-permutation null | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 0.928 | 0.257 | 0.995 | 0.696 | 0.574 | 1.5 | 0.557 | 0.659 | 0.952 | 0.919 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 14.2 | 26.8 | 8.81 | 0 | 8.52 | 23.7 | 0.013 | 0.532 | 1.34 | 0.992 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 3.27 | 3.89 | 5.11 | 9.83e-29 | 2.47 | 4.33 | 0.0939 | 0.179 | 1.99 | 1.39 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 4.14 | 5.73 | 5.01 | 1.48e-47 | 2.46 | 6.97 | 0.033 | 0.0779 | 2.29 | 1.58 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 3.8 | 3.63 | 4.99 | 9.56e-40 | 2.1 | 6.86 | 0.119 | 0.109 | 1.81 | 1.52 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 4.3 | 4.41 | 5.01 | 2.56e-51 | 2.5 | 7.4 | 0.026 | 0.0709 | 2.41 | 1.65 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 4.14 | 5.73 | 5.02 | 1.48e-47 | 2.48 | 6.9 | 0.034 | 0.0759 | 2.26 | 1.58 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | paralog_exclusion | 17458234 | 17458 | 4.09 | 5.47 | 5 | 2.39e-46 | 2.4 | 6.97 | 0.034 | 0.0789 | 2.26 | 1.58 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | class_i | 6189892 | 6.19e+03 | 2.08 | 2.67 | 2.19 | 3.68e-08 | 1.02 | 4.22 | 0.025 | 0.0529 | 2.38 | 1.76 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | NC2 | 9288636 | 9.29e+03 | 0.18 | 1.8 | 4.93 | 1 | 0.053 | 0.614 | 0.858 | 0.65 | 0.468 | 0.931 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | no_split_sibling | 6595619 | 6.6e+03 | 4.42 | 9.57 | 5.34 | 1.77e-23 | 1.26 | 15.5 | 0.023 | 0.153 | 3.35 | 1.48 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 5.2 | 5.2 | 5.08 | 2.35e-74 | 3.11 | 8.68 | 0.00799 | 0.02 | 3.07 | 2.11 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | PSEUDO_TEST |  | 17466667 | 17467 | 26.7 | 40.2 | 18.2 | 0 | 20.2 | 35.3 | 0.000999 | 0.000999 | 2.33 | 1.73 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 6.96 | 3.62 | 5.58 | 2.86e-64 | 0.792 | 61.2 | 0.027 | 0.04 | 2.95 | 1.85 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 0.822 | 0.513 | 1.06 | 0.814 | 0.273 | 2.47 | 0.439 | 0.531 | 1.07 | 0.944 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 25.9 | 40.8 | 10.1 | 0 | 7.13 | 93.9 | 0.022 | 0.0769 | 2.15 | 1.33 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 5.47 | 3.1 | 7.47 | 5.16e-42 | 1.28 | 23.4 | 0.04 | 0.136 | 5.07 | 1.78 |
| baseline | B0 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 0.616 | 0.256 | 0.95 | 0.972 | 0.35 | 1.08 | 0.95 | 0.94 | 0.609 | 0.637 |
| baseline | B1 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 0.77 | 1.03 | 0.929 | 0.873 | 0.493 | 1.2 | 0.398 | 0.639 | 1.19 | 0.822 |
| baseline | B2 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 1.23 | 4.66 | 1.36 | 0.178 | 0.0972 | 15.7 | 0.199 | 0.756 | 1.73 | 0.665 |
| baseline | B3 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 54.4 | 51.4 | 8.23 | 0 | 0.632 | 4.68e+03 | 0.034 | 0.034 | 521 | 6.05 |
| baseline | B4 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 1.7 | 0.513 | 1.28 | 0.00314 | 1.41 | 2.04 | 0.444 | 0.453 | 0.992 | 0.994 |
| baseline | B4e |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 12 | 8.89 | 13.7 | 9.89e-155 | 9.21 | 15.7 | 0.575 | 0.558 | 1 | 0.993 |
| baseline | B5 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 28.2 | 155 | 9.68 | 0 | 8.31 | 95.5 | 0.029 | 0.031 | 2.61 | 2.08 |
| baseline | B6 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 0.925 | 0.513 | 0.919 | 0.661 | 0.387 | 2.21 | 0.0889 | 0.13 | 1.34 | 1.31 |
| control | NC1 | NC1: residual-permutation null | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 0.513 | 0.256 | 0.779 | 0.993 | 0.197 | 1.33 | 0.942 | 0.876 | 0.59 | 0.679 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 14.8 | 18.1 | 8.3 | 2.94e-212 | 7.92 | 27.8 | 0.0979 | 0.201 | 1.58 | 1.2 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 6.8 | 4.14 | 6.9 | 9.46e-62 | 1.48 | 31.2 | 0.0599 | 0.145 | 3.2 | 1.55 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 7.23 | 4.14 | 5.89 | 1.56e-68 | 0.965 | 54.1 | 0.023 | 0.042 | 3.08 | 1.82 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 4.89 | 3.62 | 4.65 | 4e-34 | 0.632 | 37.8 | 0.0849 | 0.155 | 2.22 | 1.46 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 5.79 | 4.14 | 5.56 | 1.56e-46 | 0.761 | 44 | 0.0559 | 0.0789 | 2.52 | 1.66 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 6.96 | 3.62 | 5.57 | 2.86e-64 | 0.791 | 61.3 | 0.027 | 0.037 | 2.95 | 1.85 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | paralog_exclusion | 9157040 | 9.16e+03 | 6.96 | 3.62 | 5.57 | 2.73e-64 | 0.793 | 61.2 | 0.027 | 0.04 | 2.95 | 1.85 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | class_i | 693275 | 693 | 3.19 | 2.66 | 4.09 | 0.000573 | 0.0739 | 138 | 0.0919 | 0.0949 | 2.27 | 1.79 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | NC2 | 7114163 | 7.11e+03 | 0.0365 | 0.365 | 6.12 | 1 | 7.49e-08 | 17811 | 0.677 | 0.495 | 0.941 | 0.949 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | no_split_sibling | 3557652 | 3.56e+03 | 7.09 | 6.9 | 8.83 | 4.4e-26 | 1.98 | 25.4 | 0.042 | 0.004 | 2.77 | 2.09 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 6.48 | 4.14 | 5.8 | 8.36e-57 | 0.965 | 43.5 | 0.0669 | 0.147 | 2.63 | 1.54 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | PSEUDO_TEST |  | 9161007 | 9.16e+03 | 25.2 | 53.3 | 10.2 | 0 | 6.18 | 103 | 0.029 | 0.125 | 2.1 | 1.27 |
| control | PC2 |  | hct116 | PSEUDO_TEST |  | 6536088 | 6.54e+03 | 113 | 808 | 14.1 | 0 | 56 | 226 | 0.000999 | 0.000999 | 7.72 | 3.79 |

## cmh

| role | run | label | line | split | subset | CMH_OR(reg x S2 expr) | CMH_OR(own class x S2 expr) | CMH_OR(reg x HCT116 expr) [D6 diagnostic] | mean_AUROC | top10_share | OR_hub_excluded |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST |  | 0.918 | 2.16 | 0.919 | 0.506 | 0.272 | 3.67 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | PSEUDO_TEST |  | 0.977 | 1.78 | 1.01 | 0.508 | 0.149 | 1.45 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | PSEUDO_TEST |  | 2.43 | 4.52 | 2.04 | 0.64 | 0.107 | 16.7 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | PSEUDO_TEST |  | 1.27 | 4.37 | 1.23 | 0.497 | 0.42 | 5.66 |
| baseline | B0 |  | hct116 | PSEUDO_TEST |  | 1.04 | 1.06 | 1.04 | 0.499 | 0.0116 | 1.01 |
| baseline | B1 |  | hct116 | PSEUDO_TEST |  | 4.2 | 1.76 | 2.78 | 0.523 | 0.155 | 2.21 |
| baseline | B2 |  | hct116 | PSEUDO_TEST |  | 2.11 | 2.5 | 1.79 | 0.507 | 0.204 | 8.43 |
| baseline | B3 |  | hct116 | PSEUDO_TEST |  | 1.01 | 1.88 | 1.03 | 0.5 | 1 | 10.2 |
| baseline | B4 |  | hct116 | PSEUDO_TEST |  | 5.88 | 4.44 | 2.39 | 0.578 | 0.0063 | 4.89 |
| baseline | B4e |  | hct116 | PSEUDO_TEST |  | 4.46 | 2.94 | 3.97 | 0.705 | 0.0063 | 7.78 |
| baseline | B5 |  | hct116 | PSEUDO_TEST |  | 2.98 | 2.53 | 2.48 | 0.705 | 0.348 | 16.7 |
| baseline | B6 |  | hct116 | PSEUDO_TEST |  | 0.976 | 1.05 | 1.1 | 0.515 | 0.0143 | 0.882 |
| control | NC1 | NC1: residual-permutation null | hct116 | PSEUDO_TEST |  | 0.864 | 0.782 | 0.883 | 0.505 | 0.0665 | 0.692 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | PSEUDO_TEST |  | 2.27 | 2.53 | 1.98 | 0.628 | 0.063 | 10.8 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | PSEUDO_TEST |  | 0.782 | 1.7 | 0.819 | 0.509 | 0.335 | 3.99 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | PSEUDO_TEST |  | 0.919 | 2.18 | 0.919 | 0.506 | 0.272 | 3.67 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | PSEUDO_TEST |  | 0.874 | 1.81 | 0.827 | 0.505 | 0.248 | 3.04 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | PSEUDO_TEST |  | 0.934 | 2.27 | 0.939 | 0.506 | 0.286 | 3.8 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | PSEUDO_TEST |  | 0.918 | 2.17 | 0.92 | 0.506 | 0.273 | 3.74 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | paralog_exclusion | 0.912 | 2.15 | 0.913 | 0.506 | 0.272 | 3.61 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | class_i | 1.04 | 2.24 | 1.08 | 0.509 | 0.312 | 1.87 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | NC2 | 0 | 0 | 0 | 0.501 | 1 | 0.382 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | no_split_sibling | 0.866 | 3.12 | 0.857 | 0.503 | 0.414 | 3.59 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | PSEUDO_TEST |  | 1.2 | 2.77 | 1.22 | 0.501 | 0.271 | 5.34 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | PSEUDO_TEST |  | 2.55 | 4.63 | 2.12 | 0.639 | 0.11 | 16.7 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST |  | 0.987 | 3.44 | 1.03 | 0.499 | 0.891 | 3.15 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | PSEUDO_TEST |  | 0.877 | 0.908 | 0.918 | 0.491 | 0.251 | 0.502 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | PSEUDO_TEST |  | 2.04 | 4.81 | 2.55 | 0.528 | 0.632 | 20.3 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | PSEUDO_TEST |  | 0.916 | 3.52 | 0.953 | 0.498 | 0.983 | 1.51 |
| baseline | B0 |  | hek293t | PSEUDO_TEST |  | 0.635 | 0.612 | 0.629 | 0.501 | 0.022 | 0.753 |
| baseline | B1 |  | hek293t | PSEUDO_TEST |  | 0.987 | 0.909 | 0.903 | 0.495 | 0.281 | 1.76 |
| baseline | B2 |  | hek293t | PSEUDO_TEST |  | 1.45 | 0.988 | 1.41 | 0.506 | 0.345 | 2.9 |
| baseline | B3 |  | hek293t | PSEUDO_TEST |  | 0.963 | 26.9 | 0.989 | 0.5 | 1 | 47.9 |
| baseline | B4 |  | hek293t | PSEUDO_TEST |  | 1.12 | 1.11 | 0.842 | 0.479 | 0.012 | 2.77 |
| baseline | B4e |  | hek293t | PSEUDO_TEST |  | 5.71 | 3.41 | 6.08 | 0.66 | 0.012 | 10.5 |
| baseline | B5 |  | hek293t | PSEUDO_TEST |  | 1.77 | 5.64 | 2.19 | 0.66 | 0.981 | 26.7 |
| baseline | B6 |  | hek293t | PSEUDO_TEST |  | 0.943 | 1 | 0.927 | 0.501 | 0.0301 | 1.26 |
| control | NC1 | NC1: residual-permutation null | hek293t | PSEUDO_TEST |  | 0.55 | 0.456 | 0.622 | 0.506 | 0.0841 | 0.628 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | PSEUDO_TEST |  | 1.55 | 2.56 | 1.79 | 0.525 | 0.341 | 8.02 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | PSEUDO_TEST |  | 0.853 | 3.65 | 0.864 | 0.501 | 0.981 | 3.79 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | PSEUDO_TEST |  | 0.976 | 3.62 | 1.01 | 0.499 | 0.925 | 3.66 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | PSEUDO_TEST |  | 0.862 | 2.49 | 0.905 | 0.499 | 0.816 | 3.03 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | PSEUDO_TEST |  | 0.864 | 2.88 | 0.896 | 0.499 | 0.897 | 3.15 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | PSEUDO_TEST |  | 0.989 | 3.44 | 1.03 | 0.499 | 0.892 | 3.15 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | paralog_exclusion | 0.988 | 3.44 | 1.03 | 0.499 | 0.892 | 3.15 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | class_i | 0.772 | 3.28 | 0.818 | 0.5 | 0.967 | 5.21 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | NC2 | 0 | 0 | 0 | 0.499 | 1 | 0.0966 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | no_split_sibling | 0.905 | 3.17 | 0.942 | 0.5 | 0.988 | 2.82 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | PSEUDO_TEST |  | 0.935 | 3.18 | 0.946 | 0.5 | 0.893 | 3.28 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | PSEUDO_TEST |  | 2.04 | 4.65 | 2.55 | 0.528 | 0.63 | 16.7 |
| control | PC2 |  | hct116 | PSEUDO_TEST |  | 13.2 | 30.7 | 12.4 | 0.569 | 0.625 | 168 |

## sign

| role | run | label | line | split | subset | sign_concordant_OR | n_overlap | concordance | c0 | sign_OR | sign_perm_p | reg_weighted_concordance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST |  | 5.32 | 158 | 0.551 | 0.467 | 2.25 | 0.0879 | 0.526 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | PSEUDO_TEST |  | 1.95 | 66 | 0.545 | 0.484 | 1.71 | 0.493 | 0.532 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | PSEUDO_TEST |  | 26.3 | 940 | 0.44 | 0.413 | 1.39 | 0.000999 | 0.471 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | PSEUDO_TEST |  | 6.37 | 266 | 0.35 | 0.363 | 0.815 | 0.121 | 0.461 |
| baseline | B0 |  | hct116 | PSEUDO_TEST |  | 1.07 | 41 | 0.61 | 0.61 | 1.55 | 1 |  |
| baseline | B1 |  | hct116 | PSEUDO_TEST |  | 1.27 | 49 | 0.51 | 0.6 | 0.106 | 1 | 0.771 |
| baseline | B2 |  | hct116 | PSEUDO_TEST |  | 7.52 | 206 | 0.825 | 0.825 | 4.67 | 1 | 0.728 |
| baseline | B3 |  | hct116 | PSEUDO_TEST |  | 6.36 | 178 | 0.826 | 0.826 | 4.68 | 1 | 0.826 |
| baseline | B4 |  | hct116 | PSEUDO_TEST |  | 3.26 | 207 | 0.367 | 0.367 | 0.582 | 1 | 0.25 |
| baseline | B4e |  | hct116 | PSEUDO_TEST |  | 7.89 | 317 | 0.574 | 0.574 | 1.35 | 1 | 0.536 |
| baseline | B5 |  | hct116 | PSEUDO_TEST |  | 18.4 | 608 | 0.683 | 0.683 | 2.15 | 1 | 0.612 |
| baseline | B6 |  | hct116 | PSEUDO_TEST |  | 0.906 | 35 | 0.486 | 0.469 | 1.18 | 0.327 | 0.4 |
| control | NC1 | NC1: residual-permutation null | hct116 | PSEUDO_TEST |  | 1.07 | 36 | 0.583 | 0.543 | 1.48 | 1 | 0.75 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | PSEUDO_TEST |  | 13.5 | 529 | 0.437 | 0.486 | 0.617 | 0.564 | 0.521 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | PSEUDO_TEST |  | 3.85 | 126 | 0.516 | 0.461 | 1.79 | 0.782 | 0.381 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | PSEUDO_TEST |  | 5.38 | 159 | 0.553 | 0.467 | 2.29 | 0.0709 | 0.526 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | PSEUDO_TEST |  | 4.91 | 146 | 0.555 | 0.459 | 2.73 | 0.004 | 0.581 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | PSEUDO_TEST |  | 5.46 | 165 | 0.539 | 0.461 | 2.16 | 0.143 | 0.46 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | PSEUDO_TEST |  | 5.38 | 159 | 0.553 | 0.467 | 2.3 | 0.0809 | 0.526 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | paralog_exclusion | 5.26 | 157 | 0.548 | 0.466 | 2.2 | 0.0859 | 0.526 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | class_i | 2.26 | 70 | 0.5 | 0.473 | 1.28 | 0.859 | 0.48 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | NC2 | 3.89 | 0 |  |  |  |  |  |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | no_split_sibling | 6.36 | 70 | 0.643 | 0.493 | 3.92 | 0.00599 | 0.575 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | PSEUDO_TEST |  | 5.04 | 199 | 0.412 | 0.451 | 0.685 | 0.885 | 0.41 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | PSEUDO_TEST |  | 26.3 | 958 | 0.43 | 0.436 | 0.939 | 0.632 | 0.437 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST |  | 23.3 | 133 | 0.481 | 0.49 | 0.907 | 0.966 | 0.363 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | PSEUDO_TEST |  | 0.707 | 16 | 0.438 | 0.453 | 0.75 | 0.508 | 0.4 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | PSEUDO_TEST |  | 122 | 468 | 0.637 | 0.499 | 3.14 | 0.004 | 0.609 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | PSEUDO_TEST |  | 27.1 | 105 | 0.457 | 0.443 | 1.19 | 0.786 | 0.415 |
| baseline | B0 |  | hek293t | PSEUDO_TEST |  | 0.75 | 12 | 0.583 | 0.583 | 1.36 | 1 |  |
| baseline | B1 |  | hek293t | PSEUDO_TEST |  | 1.25 | 15 | 0.8 | 0.533 | 16 | 0.396 | 0.8 |
| baseline | B2 |  | hek293t | PSEUDO_TEST |  | 1.97 | 24 | 0.833 | 0.833 | 4.56 | 1 | 1 |
| baseline | B3 |  | hek293t | PSEUDO_TEST |  | 85.7 | 911 | 0.753 | 0.753 | 3.04 | 1 | 0.753 |
| baseline | B4 |  | hek293t | PSEUDO_TEST |  | 2.8 | 33 | 0.788 | 0.788 | 3.53 | 1 |  |
| baseline | B4e |  | hek293t | PSEUDO_TEST |  | 14 | 226 | 0.558 | 0.558 | 1.26 | 1 | 0.54 |
| baseline | B5 |  | hek293t | PSEUDO_TEST |  | 29.8 | 506 | 0.52 | 0.52 | 1.08 | 1 | 0.664 |
| baseline | B6 |  | hek293t | PSEUDO_TEST |  | 1 | 18 | 0.556 | 0.481 | 2 | 1 |  |
| control | NC1 | NC1: residual-permutation null | hek293t | PSEUDO_TEST |  | 0.401 | 10 | 0.4 | 0.5 | 0.375 | 1 |  |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | PSEUDO_TEST |  | 69.5 | 277 | 0.617 | 0.518 | 2.35 | 1 | 0.617 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | PSEUDO_TEST |  | 40 | 130 | 0.746 | 0.557 | 9 | 0.000999 | 0.723 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | PSEUDO_TEST |  | 24.8 | 138 | 0.493 | 0.492 | 1 | 0.933 | 0.392 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | PSEUDO_TEST |  | 19.7 | 94 | 0.596 | 0.488 | 3.88 | 0.227 | 0.424 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | PSEUDO_TEST |  | 22.2 | 111 | 0.532 | 0.479 | 2.17 | 0.467 | 0.4 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | PSEUDO_TEST |  | 23.4 | 133 | 0.481 | 0.49 | 0.907 | 0.966 | 0.363 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | paralog_exclusion | 23.3 | 133 | 0.481 | 0.49 | 0.907 | 0.966 | 0.363 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | class_i | 2.5 | 12 | 0.417 | 0.333 | 2.33 | 0.413 | 0.433 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | NC2 | 0.632 | 0 |  |  |  |  |  |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | no_split_sibling | 18.9 | 51 | 0.49 | 0.461 | 1.43 | 0.431 | 0.388 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | PSEUDO_TEST |  | 25.2 | 124 | 0.54 | 0.443 | 3.15 | 0.493 | 0.559 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | PSEUDO_TEST |  | 93.7 | 457 | 0.503 | 0.5 | 1.03 | 0.766 | 0.452 |
| control | PC2 |  | hct116 | PSEUDO_TEST |  | 4.78e+03 | 374 | 0.979 | 0.522 | 2.13e+03 | 0.000999 | 0.981 |

## callset

| role | run | label | line | split | subset | n_calls | callset_OR | callset_fisher_p |
|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST |  | 997722 | 1.23 | 9.22e-24 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | PSEUDO_TEST |  | 997722 | 1.23 | 9.22e-24 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | PSEUDO_TEST |  | 997722 | 1.23 | 9.22e-24 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | PSEUDO_TEST |  | 0 | 449 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | PSEUDO_TEST |  | 279030 | 1.06 | 0.0599 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | PSEUDO_TEST |  | 997292 | 1.23 | 3.87e-23 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | PSEUDO_TEST |  | 2785276 | 1.24 | 1.82e-60 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | PSEUDO_TEST |  | 923243 | 1.23 | 1.26e-22 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | PSEUDO_TEST |  | 991048 | 1.23 | 8.75e-24 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | paralog_exclusion | 997231 | 1.23 | 9.04e-24 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | class_i | 357729 | 1.14 | 9.06e-09 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | NC2 | 525971 | 0.962 | 0.691 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | no_split_sibling | 362403 | 1.56 | 5.79e-48 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST |  | 494245 | 0.849 | 1 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | PSEUDO_TEST |  | 494245 | 0.849 | 1 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | PSEUDO_TEST |  | 494245 | 0.849 | 1 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | PSEUDO_TEST |  | 0 | 470 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | PSEUDO_TEST |  | 140579 | 0.862 | 0.993 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | PSEUDO_TEST |  | 494000 | 0.846 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | PSEUDO_TEST |  | 1384219 | 0.826 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | PSEUDO_TEST |  | 456304 | 0.813 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | PSEUDO_TEST |  | 490960 | 0.847 | 1 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | paralog_exclusion | 494005 | 0.849 | 1 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | class_i | 31490 | 0.915 | 0.873 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | NC2 | 383792 | 0.786 | 1 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | no_split_sibling | 181603 | 0.826 | 1 |
| control | PC2 |  | hct116 | PSEUDO_TEST |  | 12839 | 70.9 | 0 |

## C1

| role | run | label | line | split | subset | jk_CI95_lo | jk_CI95_hi | dlogOR | OR_run | OR_C1 | dOR | c1_mode |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| control | g1_t100_w0__a3_r095_s2_w0 | C1 pair-specific gain | hct116 | PSEUDO_TEST |  | -0.535 | 0.0669 | -0.234 | 4.11 | 5.2 | -1.08 | score (approximation) |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | C1 pair-specific gain | hct116 | PSEUDO_TEST |  | -0.168 | 0.127 | -0.0205 | 26.2 | 26.7 | -0.543 | pair_score x target_weight |
| control | g1_t100_w0__a3_r095_s2_w0 | C1 pair-specific gain | hek293t | PSEUDO_TEST |  | -0.223 | 0.366 | 0.0715 | 6.96 | 6.48 | 0.48 | score (approximation) |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | C1 pair-specific gain | hek293t | PSEUDO_TEST |  | -0.129 | 0.18 | 0.0256 | 25.9 | 25.2 | 0.655 | pair_score x target_weight |

## paired

| role | run | label | line | split | subset | jk_CI95_lo | jk_CI95_hi | dlogOR | boot_lb90 | dlogCMH | boot_cmh_ub90 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| comparison | a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs a3_r095_s2_w0 | hct116 | PSEUDO_TEST |  | 0.501 | 1.26 | 0.881 | 0.589 | -0.0622 | 0.282 |
| comparison | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | hct116 | PSEUDO_TEST |  | -2.21 | -1.49 | -1.85 | -2.23 | -0.973 | -0.687 |
| comparison | focal__g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs focal__g1_t100_w0__a3_r095_s2_w0 | hct116 | PSEUDO_TEST |  | -2.03 | 0.972 | -0.53 | -1.52 | -0.321 | 0.136 |
| comparison | B0 | g1_t100_w0__a3_r095_s2_w0 vs B0 | hct116 | PSEUDO_TEST |  | 0.735 | 1.98 | 1.36 | 0.727 | -0.122 | 0.501 |
| comparison | B1 | g1_t100_w0__a3_r095_s2_w0 vs B1 | hct116 | PSEUDO_TEST |  | 0.265 | 2.09 | 1.18 | 0.333 | -1.52 | -0.893 |
| comparison | B2 | g1_t100_w0__a3_r095_s2_w0 vs B2 | hct116 | PSEUDO_TEST |  | -0.856 | 0.318 | -0.269 | -0.783 | -0.832 | -0.318 |
| comparison | B3 | g1_t100_w0__a3_r095_s2_w0 vs B3 | hct116 | PSEUDO_TEST |  | -9.29 | 9.05 | -0.121 | -5.02 | -0.0981 |  |
| comparison | B4 | g1_t100_w0__a3_r095_s2_w0 vs B4 | hct116 | PSEUDO_TEST |  | -0.805 | 0.257 | -0.274 | -0.983 | -1.86 |  |
| comparison | B4e | g1_t100_w0__a3_r095_s2_w0 vs B4e | hct116 | PSEUDO_TEST |  | -1.16 | -0.256 | -0.71 | -1.23 | -1.58 | -1.18 |
| comparison | B5 | g1_t100_w0__a3_r095_s2_w0 vs B5 | hct116 | PSEUDO_TEST |  | -2.52 | -0.252 | -1.39 | -2.26 | -1.18 |  |
| comparison | B6 | g1_t100_w0__a3_r095_s2_w0 vs B6 | hct116 | PSEUDO_TEST |  | 0.916 | 2.12 | 1.52 | 0.904 | -0.0615 | 0.815 |
| comparison | NC1 | g1_t100_w0__a3_r095_s2_w0 vs NC1 | hct116 | PSEUDO_TEST |  | 0.892 | 2.09 | 1.49 | 0.833 | 0.0608 | 0.765 |
| comparison | C2_g2 | g1_t100_w0__a3_r095_s2_w0 vs C2_g2 | hct116 | PSEUDO_TEST |  | -1.85 | -0.626 | -1.24 | -1.83 | -0.905 | -0.508 |
| comparison | g1_t100_w0__a3_r095_s2_w0__S1 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__S1 | hct116 | PSEUDO_TEST |  | -0.443 | 0.901 | 0.229 | -0.326 | 0.16 | 0.419 |
| comparison | g1_t100_w0__a3_r095_s2_w0__noHCT116 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__noHCT116 | hct116 | PSEUDO_TEST |  | -0.022 | 0.00923 | -0.00639 | -0.0298 | -0.00146 | 0.019 |
| comparison | g1_t100_w0__a3_r095_s2_w0__hidden_factors | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__hidden_factors | hct116 | PSEUDO_TEST |  | -0.226 | 0.386 | 0.08 | -0.187 | 0.0491 | 0.376 |
| comparison | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion | hct116 | PSEUDO_TEST |  | -0.118 | 0.0303 | -0.0439 | -0.142 | -0.0176 | 0.08 |
| comparison | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop | hct116 | PSEUDO_TEST |  | -0.0302 | 0.0174 | -0.00639 | -0.0464 | -0.000147 | 0.0136 |
| comparison | a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs a3_r095_s2_w0 | hek293t | PSEUDO_TEST |  | 0.753 | 3.52 | 2.14 | -0.909 | 0.118 |  |
| comparison | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | hek293t | PSEUDO_TEST |  | -2.44 | -0.189 | -1.31 | -4 | -0.726 |  |
| comparison | focal__g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs focal__g1_t100_w0__a3_r095_s2_w0 | hek293t | PSEUDO_TEST |  | -0.634 | 1.12 | 0.241 | -1.39 | 0.0743 |  |
| comparison | B0 | g1_t100_w0__a3_r095_s2_w0 vs B0 | hek293t | PSEUDO_TEST |  | 0.339 | 4.51 | 2.43 | -1.12 | 0.44 |  |
| comparison | B1 | g1_t100_w0__a3_r095_s2_w0 vs B1 | hek293t | PSEUDO_TEST |  | 0.169 | 4.23 | 2.2 | -1.61 | -0.000437 |  |
| comparison | B2 | g1_t100_w0__a3_r095_s2_w0 vs B2 | hek293t | PSEUDO_TEST |  | -1.7 | 5.16 | 1.73 | -0.748 | -0.387 |  |
| comparison | B3 | g1_t100_w0__a3_r095_s2_w0 vs B3 | hek293t | PSEUDO_TEST |  | -6.58 | 2.47 | -2.06 | -6.48 | 0.0245 |  |
| comparison | B4 | g1_t100_w0__a3_r095_s2_w0 vs B4 | hek293t | PSEUDO_TEST |  | -0.798 | 3.62 | 1.41 | -2.3 | -0.129 |  |
| comparison | B4e | g1_t100_w0__a3_r095_s2_w0 vs B4e | hek293t | PSEUDO_TEST |  | -2.67 | 1.58 | -0.545 | -4.4 | -1.76 |  |
| comparison | B5 | g1_t100_w0__a3_r095_s2_w0 vs B5 | hek293t | PSEUDO_TEST |  | -3.62 | 0.824 | -1.4 | -5.35 | -0.582 |  |
| comparison | B6 | g1_t100_w0__a3_r095_s2_w0 vs B6 | hek293t | PSEUDO_TEST |  | -0.34 | 4.38 | 2.02 | -1.34 | 0.0449 |  |
| comparison | NC1 | g1_t100_w0__a3_r095_s2_w0 vs NC1 | hek293t | PSEUDO_TEST |  | -0.0653 | 5.28 | 2.61 | -1.42 | 0.584 |  |
| comparison | C2_g2 | g1_t100_w0__a3_r095_s2_w0 vs C2_g2 | hek293t | PSEUDO_TEST |  | -2.42 | 0.909 | -0.757 | -4.47 | -0.451 |  |
| comparison | g1_t100_w0__a3_r095_s2_w0__S1 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__S1 | hek293t | PSEUDO_TEST |  | -0.678 | 0.725 | 0.0233 | -0.547 | 0.146 |  |
| comparison | g1_t100_w0__a3_r095_s2_w0__noHCT116 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__noHCT116 | hek293t | PSEUDO_TEST |  | -0.244 | 0.168 | -0.0377 | -0.523 | 0.0105 |  |
| comparison | g1_t100_w0__a3_r095_s2_w0__hidden_factors | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__hidden_factors | hek293t | PSEUDO_TEST |  | 0.0374 | 0.669 | 0.353 | -0.491 | 0.135 |  |
| comparison | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion | hek293t | PSEUDO_TEST |  | 0.0305 | 0.338 | 0.184 | -0.0228 | 0.133 |  |
| comparison | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop | hek293t | PSEUDO_TEST |  | -0.0029 | 0.0029 | 0 | -0.00223 | -0.00228 |  |

## runtime

| role | run | label | line | split | subset | eval_wall_s | score_wall_s | threads |
|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST |  | 2.49 | 67.2 | {'cpu_count': 10, 'blas': []} |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | PSEUDO_TEST |  | 2.68 | 66.1 | {'cpu_count': 10, 'blas': []} |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | PSEUDO_TEST |  | 3.23 | 82.7 | {'cpu_count': 10, 'blas': []} |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | PSEUDO_TEST |  | 2.51 | 51.5 | {'cpu_count': 10, 'blas': []} |
| baseline | B0 |  | hct116 | PSEUDO_TEST |  | 2.09 | 2.73 | {'cpu_count': 10, 'blas': []} |
| baseline | B1 |  | hct116 | PSEUDO_TEST |  | 2.12 | 4.12 | {'cpu_count': 10, 'blas': []} |
| baseline | B2 |  | hct116 | PSEUDO_TEST |  | 2.09 | 4.8 | {'cpu_count': 10, 'blas': []} |
| baseline | B3 |  | hct116 | PSEUDO_TEST |  | 1.72 | 2.39 | {'cpu_count': 10, 'blas': []} |
| baseline | B4 |  | hct116 | PSEUDO_TEST |  | 2.09 | 2.68 | {'cpu_count': 10, 'blas': []} |
| baseline | B4e |  | hct116 | PSEUDO_TEST |  | 2.23 | 2.86 | {'cpu_count': 10, 'blas': []} |
| baseline | B5 |  | hct116 | PSEUDO_TEST |  | 1.91 | 2.41 | {'cpu_count': 10, 'blas': []} |
| baseline | B6 |  | hct116 | PSEUDO_TEST |  | 2.22 | 12.2 | {'cpu_count': 10, 'blas': []} |
| control | NC1 | NC1: residual-permutation null | hct116 | PSEUDO_TEST |  | 2.33 | 27.5 | {'cpu_count': 10, 'blas': []} |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | PSEUDO_TEST |  | 1.91 | 42.8 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | PSEUDO_TEST |  | 2.56 | 53.9 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | PSEUDO_TEST |  | 2.68 | 66.5 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | PSEUDO_TEST |  | 2.71 | 90.7 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | PSEUDO_TEST |  | 2.67 | 66.6 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | PSEUDO_TEST |  | 2.64 | 66.1 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | paralog_exclusion | 2.74 | 67.2 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | class_i | 1.9 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | NC2 | 1.9 | 67.2 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | no_split_sibling | 1.89 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | PSEUDO_TEST |  | 0.929 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | PSEUDO_TEST |  | 1 | 82.7 | {'cpu_count': 10, 'blas': []} |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST |  | 1.84 | 67.2 | {'cpu_count': 10, 'blas': []} |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | PSEUDO_TEST |  | 2.15 | 66.1 | {'cpu_count': 10, 'blas': []} |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | PSEUDO_TEST |  | 2.35 | 82.7 | {'cpu_count': 10, 'blas': []} |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | PSEUDO_TEST |  | 2.07 | 51.5 | {'cpu_count': 10, 'blas': []} |
| baseline | B0 |  | hek293t | PSEUDO_TEST |  | 1.62 | 2.73 | {'cpu_count': 10, 'blas': []} |
| baseline | B1 |  | hek293t | PSEUDO_TEST |  | 1.61 | 4.12 | {'cpu_count': 10, 'blas': []} |
| baseline | B2 |  | hek293t | PSEUDO_TEST |  | 1.68 | 4.8 | {'cpu_count': 10, 'blas': []} |
| baseline | B3 |  | hek293t | PSEUDO_TEST |  | 1.48 | 2.39 | {'cpu_count': 10, 'blas': []} |
| baseline | B4 |  | hek293t | PSEUDO_TEST |  | 1.68 | 2.68 | {'cpu_count': 10, 'blas': []} |
| baseline | B4e |  | hek293t | PSEUDO_TEST |  | 1.67 | 2.86 | {'cpu_count': 10, 'blas': []} |
| baseline | B5 |  | hek293t | PSEUDO_TEST |  | 1.62 | 2.41 | {'cpu_count': 10, 'blas': []} |
| baseline | B6 |  | hek293t | PSEUDO_TEST |  | 1.77 | 12.2 | {'cpu_count': 10, 'blas': []} |
| control | NC1 | NC1: residual-permutation null | hek293t | PSEUDO_TEST |  | 1.65 | 27.5 | {'cpu_count': 10, 'blas': []} |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | PSEUDO_TEST |  | 1.52 | 42.8 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | PSEUDO_TEST |  | 2.07 | 53.9 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | PSEUDO_TEST |  | 2.09 | 66.5 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | PSEUDO_TEST |  | 2.12 | 90.7 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | PSEUDO_TEST |  | 2.08 | 66.6 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | PSEUDO_TEST |  | 2.1 | 66.1 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | paralog_exclusion | 2.15 | 67.2 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | class_i | 1.47 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | NC2 | 1.73 | 67.2 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | no_split_sibling | 1.61 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | PSEUDO_TEST |  | 0.515 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | PSEUDO_TEST |  | 0.543 | 82.7 | {'cpu_count': 10, 'blas': []} |
| control | PC2 |  | hct116 | PSEUDO_TEST |  | 1.8 |  |  |

## class:RBP

| role | run | label | line | split | subset | OR@0.1% | fisher_p(indep.) | n_regulators | selected | positives | a |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST |  | 3.67 | 3.01e-06 | 13 | 149 | 5.13e+03 | 20 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | PSEUDO_TEST |  | 4.83 | 0.000265 | 13 | 53 | 5.13e+03 | 9 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | PSEUDO_TEST |  | 13.5 | 1.49e-78 | 13 | 334 | 5.13e+03 | 120 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | PSEUDO_TEST |  | 1.23 | 0.35 | 13 | 141 | 5.13e+03 | 7 |
| baseline | B0 |  | hct116 | PSEUDO_TEST |  | 1.59 | 0.13 | 13 | 143 | 5.13e+03 | 9 |
| baseline | B1 |  | hct116 | PSEUDO_TEST |  | 0.0343 | 1 | 13 | 342 | 5.13e+03 | 0 |
| baseline | B2 |  | hct116 | PSEUDO_TEST |  | 1.07 | 0.386 | 13 | 644 | 5.13e+03 | 28 |
| baseline | B3 |  | hct116 | PSEUDO_TEST |  | 23.6 | 1 | 13 | 0 | 5.13e+03 | 0 |
| baseline | B4 |  | hct116 | PSEUDO_TEST |  | 2.88 | 0.000809 | 13 | 129 | 5.13e+03 | 14 |
| baseline | B4e |  | hct116 | PSEUDO_TEST |  | 5.8 | 5.76e-11 | 13 | 127 | 5.13e+03 | 25 |
| baseline | B5 |  | hct116 | PSEUDO_TEST |  | 0.0223 | 1 | 13 | 526 | 5.13e+03 | 0 |
| baseline | B6 |  | hct116 | PSEUDO_TEST |  | 0.761 | 0.769 | 13 | 128 | 5.13e+03 | 4 |
| control | NC1 | NC1: residual-permutation null | hct116 | PSEUDO_TEST |  | 1.43 | 0.211 | 13 | 140 | 5.13e+03 | 8 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | PSEUDO_TEST |  | 4.45 | 6.75e-14 | 13 | 266 | 5.13e+03 | 42 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | PSEUDO_TEST |  | 3.82 | 9.57e-08 | 13 | 180 | 5.13e+03 | 25 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | PSEUDO_TEST |  | 3.64 | 3.34e-06 | 13 | 150 | 5.13e+03 | 20 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | PSEUDO_TEST |  | 1.25 | 0.23 | 13 | 319 | 5.13e+03 | 16 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | PSEUDO_TEST |  | 3.75 | 2.19e-06 | 13 | 146 | 5.13e+03 | 20 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | PSEUDO_TEST |  | 3.64 | 3.34e-06 | 13 | 150 | 5.13e+03 | 20 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | paralog_exclusion | 3.67 | 2.96e-06 | 13 | 149 | 5.12e+03 | 20 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | class_i | 5.84 | 2.52e-06 | 11 | 57 | 5.13e+03 | 13 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | NC2 | 13256 | 1 | 2 | 0 | 1 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | no_split_sibling | 4.51 | 0.00823 | 5 | 29 | 2.13e+03 | 5 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | PSEUDO_TEST |  | 2.57 | 0.00153 | 13 | 153 | 5.13e+03 | 15 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | PSEUDO_TEST |  | 12.7 | 2.48e-76 | 13 | 347 | 5.13e+03 | 120 |

## class:TF

| role | run | label | line | split | subset | OR@0.1% | fisher_p(indep.) | n_regulators | selected | positives | a |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST |  | 10.7 | 1.42e-08 | 206 | 1.79e+03 | 1.15e+03 | 11 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | PSEUDO_TEST |  | 4.03 | 0.009 | 206 | 2.14e+03 | 1.15e+03 | 5 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | PSEUDO_TEST |  | 20.1 | 1.15e-17 | 206 | 1.58e+03 | 1.15e+03 | 18 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | PSEUDO_TEST |  | 14.8 | 4.41e-26 | 206 | 3.85e+03 | 1.15e+03 | 32 |
| baseline | B0 |  | hct116 | PSEUDO_TEST |  | 0.879 | 0.68 | 206 | 1.95e+03 | 1.15e+03 | 1 |
| baseline | B1 |  | hct116 | PSEUDO_TEST |  | 1.63 | 0.459 | 206 | 1.06e+03 | 1.15e+03 | 1 |
| baseline | B2 |  | hct116 | PSEUDO_TEST |  | 1.82 | 0.423 | 206 | 943 | 1.15e+03 | 1 |
| baseline | B3 |  | hct116 | PSEUDO_TEST |  | 1.71e+03 | 1 | 206 | 0 | 1.15e+03 | 0 |
| baseline | B4 |  | hct116 | PSEUDO_TEST |  | 19.2 | 4.54e-19 | 206 | 1.84e+03 | 1.15e+03 | 20 |
| baseline | B4e |  | hct116 | PSEUDO_TEST |  | 4.24 | 0.00734 | 206 | 2.04e+03 | 1.15e+03 | 5 |
| baseline | B5 |  | hct116 | PSEUDO_TEST |  | 17.7 | 1 | 206 | 48 | 1.15e+03 | 0 |
| baseline | B6 |  | hct116 | PSEUDO_TEST |  | 0.351 | 1 | 206 | 2.44e+03 | 1.15e+03 | 0 |
| control | NC1 | NC1: residual-permutation null | hct116 | PSEUDO_TEST |  | 4.57 | 0.0293 | 206 | 1.13e+03 | 1.15e+03 | 3 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | PSEUDO_TEST |  | 8.53 | 6.77e-06 | 206 | 1.63e+03 | 1.15e+03 | 8 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | PSEUDO_TEST |  | 16.2 | 2.84e-15 | 206 | 1.84e+03 | 1.15e+03 | 17 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | PSEUDO_TEST |  | 10.7 | 1.39e-08 | 206 | 1.78e+03 | 1.15e+03 | 11 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | PSEUDO_TEST |  | 7.92 | 0.000142 | 206 | 1.31e+03 | 1.15e+03 | 6 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | PSEUDO_TEST |  | 10.8 | 2.88e-09 | 206 | 1.93e+03 | 1.15e+03 | 12 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | PSEUDO_TEST |  | 10.7 | 1.37e-08 | 206 | 1.78e+03 | 1.15e+03 | 11 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | paralog_exclusion | 10.7 | 1.42e-08 | 206 | 1.79e+03 | 1.15e+03 | 11 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | class_i | 4.72 | 0.192 | 39 | 154 | 518 | 1 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | NC2 | 0.443 | 1 | 110 | 3.47e+03 | 341 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | no_split_sibling | 11.3 | 1.93e-07 | 45 | 856 | 408 | 9 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | PSEUDO_TEST |  | 19.1 | 3.77e-18 | 206 | 1.75e+03 | 1.15e+03 | 19 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | PSEUDO_TEST |  | 11.8 | 5.3e-09 | 206 | 1.62e+03 | 1.15e+03 | 11 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST |  | 32.8 | 1 | 60 | 2 | 3.78e+03 | 0 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | PSEUDO_TEST |  | 0.249 | 1 | 60 | 328 | 3.78e+03 | 0 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | PSEUDO_TEST |  | 2.18 | 1 | 60 | 37 | 3.78e+03 | 0 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | PSEUDO_TEST |  | 18.2 | 1 | 60 | 4 | 3.78e+03 | 0 |
| baseline | B0 |  | hek293t | PSEUDO_TEST |  | 0.128 | 1 | 60 | 637 | 3.78e+03 | 0 |
| baseline | B1 |  | hek293t | PSEUDO_TEST |  | 0.207 | 1 | 60 | 395 | 3.78e+03 | 0 |
| baseline | B2 |  | hek293t | PSEUDO_TEST |  | 1.27 | 1 | 60 | 64 | 3.78e+03 | 0 |
| baseline | B3 |  | hek293t | PSEUDO_TEST |  | 164 | 1 | 60 | 0 | 3.78e+03 | 0 |
| baseline | B4 |  | hek293t | PSEUDO_TEST |  | 1.61 | 0.174 | 60 | 616 | 3.78e+03 | 6 |
| baseline | B4e |  | hek293t | PSEUDO_TEST |  | 10.1 | 4.04e-24 | 60 | 641 | 3.78e+03 | 37 |
| baseline | B5 |  | hek293t | PSEUDO_TEST |  | 164 | 1 | 60 | 0 | 3.78e+03 | 0 |
| baseline | B6 |  | hek293t | PSEUDO_TEST |  | 1.07 | 0.518 | 60 | 619 | 3.78e+03 | 4 |
| control | NC1 | NC1: residual-permutation null | hek293t | PSEUDO_TEST |  | 0.973 | 0.61 | 60 | 339 | 3.78e+03 | 2 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | PSEUDO_TEST |  | 20.2 | 4.63e-40 | 60 | 405 | 3.78e+03 | 44 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | PSEUDO_TEST |  | 164 | 1 | 60 | 0 | 3.78e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | PSEUDO_TEST |  | 32.8 | 1 | 60 | 2 | 3.78e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | PSEUDO_TEST |  | 0.385 | 1 | 60 | 212 | 3.78e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | PSEUDO_TEST |  | 23.4 | 1 | 60 | 3 | 3.78e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | PSEUDO_TEST |  | 32.8 | 1 | 60 | 2 | 3.78e+03 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | paralog_exclusion | 32.7 | 1 | 60 | 2 | 3.77e+03 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | class_i | 8.79e+03 | 1 | 3 | 0 | 3 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | NC2 | 146 | 1 | 34 | 0 | 2.39e+03 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | no_split_sibling | 982 | 1 | 17 | 0 | 180 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | PSEUDO_TEST |  | 32.8 | 1 | 60 | 2 | 3.78e+03 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | PSEUDO_TEST |  | 3.34 | 1 | 60 | 24 | 3.78e+03 | 0 |
| control | PC2 |  | hct116 | PSEUDO_TEST |  | 193 | 2.41e-40 | 46 | 1.16e+03 | 66 | 22 |

## class:kinase

| role | run | label | line | split | subset | OR@0.1% | fisher_p(indep.) | n_regulators | selected | positives | a |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST |  | 3.88 | 0.229 | 39 | 139 | 702 | 1 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | PSEUDO_TEST |  | 4.73 | 0.192 | 39 | 114 | 702 | 1 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | PSEUDO_TEST |  | 62.7 | 2.6e-47 | 39 | 338 | 702 | 34 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | PSEUDO_TEST |  | 8.75 | 1 | 39 | 30 | 702 | 0 |
| baseline | B0 |  | hct116 | PSEUDO_TEST |  | 1.32 | 0.533 | 39 | 407 | 702 | 1 |
| baseline | B1 |  | hct116 | PSEUDO_TEST |  | 1.08 | 1 | 39 | 247 | 702 | 0 |
| baseline | B2 |  | hct116 | PSEUDO_TEST |  | 0.257 | 1 | 39 | 1.04e+03 | 702 | 0 |
| baseline | B3 |  | hct116 | PSEUDO_TEST |  | 534 | 1 | 39 | 0 | 702 | 0 |
| baseline | B4 |  | hct116 | PSEUDO_TEST |  | 10.3 | 8.05e-06 | 39 | 373 | 702 | 7 |
| baseline | B4e |  | hct116 | PSEUDO_TEST |  | 11.7 | 7.46e-07 | 39 | 378 | 702 | 8 |
| baseline | B5 |  | hct116 | PSEUDO_TEST |  | 534 | 1 | 39 | 0 | 702 | 0 |
| baseline | B6 |  | hct116 | PSEUDO_TEST |  | 0.643 | 1 | 39 | 414 | 702 | 0 |
| control | NC1 | NC1: residual-permutation null | hct116 | PSEUDO_TEST |  | 0.884 | 0.678 | 39 | 605 | 702 | 1 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | PSEUDO_TEST |  | 15.4 | 4.97e-14 | 39 | 584 | 702 | 16 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | PSEUDO_TEST |  | 7.04 | 0.134 | 39 | 77 | 702 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | PSEUDO_TEST |  | 3.88 | 0.229 | 39 | 139 | 702 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | PSEUDO_TEST |  | 2.33 | 1 | 39 | 114 | 702 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | PSEUDO_TEST |  | 4.02 | 0.222 | 39 | 134 | 702 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | PSEUDO_TEST |  | 3.85 | 0.23 | 39 | 140 | 702 | 1 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | paralog_exclusion | 3.87 | 0.229 | 39 | 139 | 702 | 1 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | class_i | 3.97 | 0.226 | 16 | 60 | 655 | 1 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | NC2 | 4.61e+03 | 1 | 15 | 0 | 31 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | no_split_sibling | 1.52e+03 | 1 | 17 | 2 | 21 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | PSEUDO_TEST |  | 1.83 | 1 | 39 | 145 | 702 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | PSEUDO_TEST |  | 38 | 1.24e-25 | 39 | 325 | 702 | 21 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | PSEUDO_TEST |  | 8.42 | 1 | 15 | 91 | 101 | 0 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| baseline | B0 |  | hek293t | PSEUDO_TEST |  | 4.34 | 1 | 15 | 177 | 101 | 0 |
| baseline | B1 |  | hek293t | PSEUDO_TEST |  | 3.99 | 0.224 | 15 | 392 | 101 | 1 |
| baseline | B2 |  | hek293t | PSEUDO_TEST |  | 1.54 | 1 | 15 | 497 | 101 | 0 |
| baseline | B3 |  | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| baseline | B4 |  | hek293t | PSEUDO_TEST |  | 4.89 | 1 | 15 | 157 | 101 | 0 |
| baseline | B4e |  | hek293t | PSEUDO_TEST |  | 31.5 | 0.000148 | 15 | 155 | 101 | 3 |
| baseline | B5 |  | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| baseline | B6 |  | hek293t | PSEUDO_TEST |  | 8.45 | 0.113 | 15 | 186 | 101 | 1 |
| control | NC1 | NC1: residual-permutation null | hek293t | PSEUDO_TEST |  | 6.04 | 1 | 15 | 127 | 101 | 0 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | paralog_exclusion | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | NC2 | 1.29e+03 | 1 | 10 | 0 | 81 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | no_split_sibling | 2.05e+03 | 1 | 6 | 0 | 30 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | PSEUDO_TEST |  | 1.54e+03 | 1 | 15 | 0 | 101 | 0 |
| control | PC2 |  | hct116 | PSEUDO_TEST |  | 1.25e+03 | 9.89e-17 | 10 | 29 | 26 | 6 |

## class:other

| role | run | label | line | split | subset | OR@0.1% | fisher_p(indep.) | n_regulators | selected | positives | a |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST |  | 3.97 | 7.01e-34 | 1.54e+03 | 14385 | 30721 | 117 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | PSEUDO_TEST |  | 1.52 | 0.00492 | 1.54e+03 | 14696 | 30721 | 46 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | PSEUDO_TEST |  | 25.7 | 0 | 1.54e+03 | 14630 | 30721 | 724 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | PSEUDO_TEST |  | 8.88 | 2.18e-127 | 1.54e+03 | 12472 | 30721 | 224 |
| baseline | B0 |  | hct116 | PSEUDO_TEST |  | 0.979 | 0.57 | 1.54e+03 | 14810 | 30721 | 30 |
| baseline | B1 |  | hct116 | PSEUDO_TEST |  | 1.47 | 0.00716 | 1.54e+03 | 15821 | 30721 | 48 |
| baseline | B2 |  | hct116 | PSEUDO_TEST |  | 5.48 | 4.31e-65 | 1.54e+03 | 14747 | 30721 | 165 |
| baseline | B3 |  | hct116 | PSEUDO_TEST |  | 4.99 | 3.35e-64 | 1.54e+03 | 17467 | 30721 | 178 |
| baseline | B4 |  | hct116 | PSEUDO_TEST |  | 5.26 | 3.04e-61 | 1.54e+03 | 14983 | 30721 | 161 |
| baseline | B4e |  | hct116 | PSEUDO_TEST |  | 9.11 | 7.78e-157 | 1.54e+03 | 14789 | 30721 | 272 |
| baseline | B5 |  | hct116 | PSEUDO_TEST |  | 18.5 | 0 | 1.54e+03 | 16745 | 30721 | 608 |
| baseline | B6 |  | hct116 | PSEUDO_TEST |  | 1.05 | 0.426 | 1.54e+03 | 14330 | 30721 | 31 |
| control | NC1 | NC1: residual-permutation null | hct116 | PSEUDO_TEST |  | 0.751 | 0.938 | 1.54e+03 | 15429 | 30721 | 24 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | PSEUDO_TEST |  | 15.9 | 0 | 1.54e+03 | 14463 | 30721 | 455 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | PSEUDO_TEST |  | 2.61 | 1.95e-13 | 1.54e+03 | 14504 | 30721 | 78 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | PSEUDO_TEST |  | 4 | 1.72e-34 | 1.54e+03 | 14382 | 30721 | 118 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | PSEUDO_TEST |  | 3.87 | 1.77e-33 | 1.54e+03 | 14994 | 30721 | 119 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | PSEUDO_TEST |  | 4.22 | 5.66e-38 | 1.54e+03 | 14239 | 30721 | 123 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | PSEUDO_TEST |  | 4 | 1.75e-34 | 1.54e+03 | 14385 | 30721 | 118 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | paralog_exclusion | 3.93 | 2.7e-33 | 1.54e+03 | 14378 | 30716 | 116 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | class_i | 1.92 | 1.75e-05 | 564 | 5.45e+03 | 26556 | 51 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | NC2 | 0.289 | 1 | 837 | 5.82e+03 | 2.39e+03 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | no_split_sibling | 4.6 | 7.56e-20 | 620 | 5.44e+03 | 13429 | 56 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | PSEUDO_TEST |  | 5.31 | 7.06e-60 | 1.54e+03 | 14392 | 30721 | 156 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | PSEUDO_TEST |  | 27.1 | 0 | 1.54e+03 | 14583 | 30721 | 759 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST |  | 8.04 | 1.99e-71 | 807 | 9.06e+03 | 15577 | 133 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | PSEUDO_TEST |  | 0.983 | 0.561 | 807 | 8.74e+03 | 15577 | 16 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | PSEUDO_TEST |  | 31.6 | 0 | 807 | 8.59e+03 | 15577 | 465 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | PSEUDO_TEST |  | 6.3 | 2.39e-47 | 807 | 9.08e+03 | 15577 | 105 |
| baseline | B0 |  | hek293t | PSEUDO_TEST |  | 0.773 | 0.847 | 807 | 8.32e+03 | 15577 | 12 |
| baseline | B1 |  | hek293t | PSEUDO_TEST |  | 0.897 | 0.692 | 807 | 8.37e+03 | 15577 | 14 |
| baseline | B2 |  | hek293t | PSEUDO_TEST |  | 1.5 | 0.0368 | 807 | 8.6e+03 | 15577 | 24 |
| baseline | B3 |  | hek293t | PSEUDO_TEST |  | 62.8 | 0 | 807 | 9.16e+03 | 15577 | 911 |
| baseline | B4 |  | hek293t | PSEUDO_TEST |  | 1.74 | 0.00535 | 807 | 8.37e+03 | 15577 | 27 |
| baseline | B4e |  | hek293t | PSEUDO_TEST |  | 12.3 | 9.07e-129 | 807 | 8.34e+03 | 15577 | 185 |
| baseline | B5 |  | hek293t | PSEUDO_TEST |  | 36.6 | 0 | 807 | 8.09e+03 | 15577 | 502 |
| baseline | B6 |  | hek293t | PSEUDO_TEST |  | 0.837 | 0.774 | 807 | 8.34e+03 | 15577 | 13 |
| control | NC1 | NC1: residual-permutation null | hek293t | PSEUDO_TEST |  | 0.494 | 0.991 | 807 | 8.67e+03 | 15577 | 8 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | PSEUDO_TEST |  | 15.4 | 5.82e-182 | 807 | 8.39e+03 | 15577 | 232 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | PSEUDO_TEST |  | 7.8 | 2.08e-68 | 807 | 9.12e+03 | 15577 | 130 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | PSEUDO_TEST |  | 8.35 | 4.93e-76 | 807 | 9.06e+03 | 15577 | 138 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | PSEUDO_TEST |  | 5.79 | 2.56e-39 | 807 | 8.74e+03 | 15577 | 93 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | PSEUDO_TEST |  | 6.68 | 2.67e-52 | 807 | 9.07e+03 | 15577 | 111 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | PSEUDO_TEST |  | 8.04 | 2.02e-71 | 807 | 9.06e+03 | 15577 | 133 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | paralog_exclusion | 8.04 | 1.92e-71 | 807 | 9.06e+03 | 15573 | 133 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | class_i | 2.99 | 0.001 | 62 | 687 | 3.8e+03 | 12 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | NC2 | 0.0417 | 1 | 642 | 7.11e+03 | 11177 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | no_split_sibling | 6.8 | 2.6e-25 | 321 | 3.55e+03 | 7.12e+03 | 51 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | PSEUDO_TEST |  | 7.48 | 2.13e-63 | 807 | 9.06e+03 | 15577 | 124 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | PSEUDO_TEST |  | 30.9 | 0 | 807 | 8.57e+03 | 15577 | 455 |
| control | PC2 |  | hct116 | PSEUDO_TEST |  | 120 | 0 | 623 | 5.34e+03 | 3.8e+03 | 346 |

## class:ribosomal

| role | run | label | line | split | subset | OR@0.1% | fisher_p(indep.) | n_regulators | selected | positives | a |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST |  | 1.17 | 0.371 | 15 | 1.01e+03 | 1.1e+03 | 9 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | PSEUDO_TEST |  | 1.41 | 0.285 | 15 | 462 | 1.1e+03 | 5 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | PSEUDO_TEST |  | 10.8 | 4.48e-29 | 15 | 589 | 1.1e+03 | 44 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | PSEUDO_TEST |  | 0.396 | 0.981 | 15 | 976 | 1.1e+03 | 3 |
| baseline | B0 |  | hct116 | PSEUDO_TEST |  | 0.417 | 1 | 15 | 154 | 1.1e+03 | 0 |
| baseline | B1 |  | hct116 | PSEUDO_TEST |  | 25.8 | 1 | 15 | 2 | 1.1e+03 | 0 |
| baseline | B2 |  | hct116 | PSEUDO_TEST |  | 18.4 | 1.6e-11 | 15 | 97 | 1.1e+03 | 12 |
| baseline | B3 |  | hct116 | PSEUDO_TEST |  | 129 | 1 | 15 | 0 | 1.1e+03 | 0 |
| baseline | B4 |  | hct116 | PSEUDO_TEST |  | 4.62 | 0.00551 | 15 | 145 | 1.1e+03 | 5 |
| baseline | B4e |  | hct116 | PSEUDO_TEST |  | 6.98 | 0.000102 | 15 | 137 | 1.1e+03 | 7 |
| baseline | B5 |  | hct116 | PSEUDO_TEST |  | 0.433 | 1 | 15 | 148 | 1.1e+03 | 0 |
| baseline | B6 |  | hct116 | PSEUDO_TEST |  | 0.414 | 1 | 15 | 155 | 1.1e+03 | 0 |
| control | NC1 | NC1: residual-permutation null | hct116 | PSEUDO_TEST |  | 0.398 | 1 | 15 | 161 | 1.1e+03 | 0 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | PSEUDO_TEST |  | 1.99 | 0.0535 | 15 | 527 | 1.1e+03 | 8 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | PSEUDO_TEST |  | 0.748 | 0.795 | 15 | 865 | 1.1e+03 | 5 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | PSEUDO_TEST |  | 1.16 | 0.378 | 15 | 1.01e+03 | 1.1e+03 | 9 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | PSEUDO_TEST |  | 0.891 | 0.66 | 15 | 728 | 1.1e+03 | 5 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | PSEUDO_TEST |  | 1.16 | 0.38 | 15 | 1.01e+03 | 1.1e+03 | 9 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | PSEUDO_TEST |  | 1.16 | 0.376 | 15 | 1.01e+03 | 1.1e+03 | 9 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | paralog_exclusion | 1.17 | 0.371 | 15 | 1.01e+03 | 1.1e+03 | 9 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | class_i | 1.04 | 0.538 | 14 | 467 | 1.1e+03 | 4 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | NC2 | 19255 | 1 | 1 | 0 | 0 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | PSEUDO_TEST | no_split_sibling | 2.54 | 1 | 3 | 268 | 20 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | PSEUDO_TEST |  | 1.15 | 0.39 | 15 | 1.02e+03 | 1.1e+03 | 9 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | PSEUDO_TEST |  | 11.6 | 3.31e-32 | 15 | 590 | 1.1e+03 | 47 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST |  | 12.6 | 1 | 2 | 96 | 8 | 0 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | PSEUDO_TEST |  | 272 | 1 | 2 | 4 | 8 | 0 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | PSEUDO_TEST |  | 22.9 | 0.00085 | 2 | 533 | 8 | 3 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | PSEUDO_TEST |  | 15.9 | 1 | 2 | 76 | 8 | 0 |
| baseline | B0 |  | hek293t | PSEUDO_TEST |  | 49.9 | 1 | 2 | 24 | 8 | 0 |
| baseline | B1 |  | hek293t | PSEUDO_TEST |  | 2.45e+03 | 1 | 2 | 0 | 8 | 0 |
| baseline | B2 |  | hek293t | PSEUDO_TEST |  | 272 | 1 | 2 | 4 | 8 | 0 |
| baseline | B3 |  | hek293t | PSEUDO_TEST |  | 2.45e+03 | 1 | 2 | 0 | 8 | 0 |
| baseline | B4 |  | hek293t | PSEUDO_TEST |  | 59.6 | 1 | 2 | 20 | 8 | 0 |
| baseline | B4e |  | hek293t | PSEUDO_TEST |  | 141 | 0.00843 | 2 | 22 | 8 | 1 |
| baseline | B5 |  | hek293t | PSEUDO_TEST |  | 18.5 | 0.000409 | 2 | 1.07e+03 | 8 | 4 |
| baseline | B6 |  | hek293t | PSEUDO_TEST |  | 62.7 | 1 | 2 | 19 | 8 | 0 |
| control | NC1 | NC1: residual-permutation null | hek293t | PSEUDO_TEST |  | 56.9 | 1 | 2 | 21 | 8 | 0 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | PSEUDO_TEST |  | 8.02 | 0.132 | 2 | 365 | 8 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | PSEUDO_TEST |  | 29.4 | 1 | 2 | 41 | 8 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | PSEUDO_TEST |  | 12 | 1 | 2 | 101 | 8 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | PSEUDO_TEST |  | 14.1 | 0.078 | 2 | 210 | 8 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | PSEUDO_TEST |  | 13.3 | 1 | 2 | 91 | 8 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | PSEUDO_TEST |  | 12.8 | 1 | 2 | 95 | 8 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | paralog_exclusion | 12.8 | 1 | 2 | 95 | 8 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | class_i | 188 | 1 | 2 | 6 | 8 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | PSEUDO_TEST | no_split_sibling | 305 | 1 | 1 | 11 | 1 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | PSEUDO_TEST |  | 12.6 | 1 | 2 | 96 | 8 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | PSEUDO_TEST |  | 12 | 0.0184 | 2 | 564 | 8 | 2 |
| control | PC2 |  | hct116 | PSEUDO_TEST |  | 3.81e+03 | 1 | 1 | 2 | 0 | 0 |
