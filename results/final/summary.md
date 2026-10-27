# Final TEST evaluation

## pooled

| role | run | label | line | split | subset | n_pairs | K | OR@0.1% | OR@0.01% | OR@1% | fisher_p(indep.) | jk_CI95_lo | jk_CI95_hi | perm_p | perm_p_degree | OR/median_null | OR/median_null_degree |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST |  | 17692753 | 17693 | 3.85 | 2.92 | 4.77 | 2.04e-36 | 1.33 | 11.2 | 0.0679 | 0.0609 | 2.27 | 2.1 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | TEST |  | 17692753 | 17693 | 1.46 | 0.873 | 2.02 | 0.00704 | 0.911 | 2.33 | 0.229 | 0.206 | 1.54 | 1.46 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | TEST |  | 17692753 | 17693 | 18.4 | 16.5 | 13.4 | 0 | 10.9 | 30.9 | 0.015 | 0.021 | 1.79 | 1.58 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | TEST |  | 17692753 | 17693 | 4.5 | 2.92 | 5.49 | 4.12e-50 | 1.19 | 17.1 | 0.0829 | 0.018 | 2.48 | 2.52 |
| baseline | B0 |  | hct116 | TEST |  | 17692753 | 17693 | 1.17 | 0.873 | 1 | 0.187 | 0.964 | 1.41 | 0.145 | 0.171 | 1.17 | 1.16 |
| baseline | B1 |  | hct116 | TEST |  | 17692753 | 17693 | 1.14 | 0.291 | 0.878 | 0.234 | 0.313 | 4.12 | 0.723 | 0.601 | 0.812 | 0.933 |
| baseline | B2 |  | hct116 | TEST |  | 17692753 | 17693 | 6.03 | 8.28 | 5.31 | 6.59e-87 | 2.71 | 13.4 | 0.0949 | 0.332 | 3.03 | 1.39 |
| baseline | B3 |  | hct116 | TEST |  | 17692753 | 17693 | 0.029 | 0.145 | 10.6 | 1 | 8.75e-11 | 9629726 | 0.666 | 0.422 | 0.467 | 0.954 |
| baseline | B4 |  | hct116 | TEST |  | 17692753 | 17693 | 5.73 | 10.4 | 3.52 | 3.35e-79 | 4.23 | 7.76 | 0.528 | 0.535 | 1.01 | 1.01 |
| baseline | B4e |  | hct116 | TEST |  | 17692753 | 17693 | 7.39 | 10.1 | 6.79 | 4.02e-124 | 5.55 | 9.83 | 0.531 | 0.585 | 0.999 | 0.993 |
| baseline | B5 |  | hct116 | TEST |  | 17692753 | 17693 | 22.8 | 15.9 | 12.6 | 0 | 3.27 | 159 | 0.043 | 0.035 | 2.89 | 2.34 |
| baseline | B6 |  | hct116 | TEST |  | 17692753 | 17693 | 0.611 | 1.17 | 0.811 | 0.994 | 0.349 | 1.07 | 0.989 | 0.941 | 0.576 | 0.644 |
| control | NC1 | NC1: residual-permutation null | hct116 | TEST |  | 17692753 | 17693 | 1.28 | 1.46 | 1.12 | 0.0631 | 0.484 | 3.4 | 0.163 | 0.119 | 1.56 | 1.62 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | TEST |  | 17692753 | 17693 | 10.6 | 14.1 | 7.54 | 4.25e-224 | 8.57 | 13.1 | 0.148 | 0.395 | 1.19 | 1.03 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | TEST |  | 17692753 | 17693 | 3.4 | 1.46 | 4.78 | 5.48e-28 | 1.64 | 7.06 | 0.105 | 0.118 | 2.2 | 1.89 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | TEST |  | 17692753 | 17693 | 3.79 | 2.92 | 4.78 | 3.01e-35 | 1.38 | 10.4 | 0.0689 | 0.0649 | 2.22 | 2.09 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | TEST |  | 17692753 | 17693 | 2.9 | 0.582 | 4.07 | 1.85e-19 | 1.05 | 7.99 | 0.188 | 0.119 | 1.76 | 1.84 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | TEST |  | 17692753 | 17693 | 3.49 | 1.75 | 4.54 | 1.31e-29 | 1.38 | 8.85 | 0.0969 | 0.0739 | 2.05 | 2.02 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | TEST |  | 17692753 | 17693 | 3.73 | 2.92 | 4.73 | 4.31e-34 | 1.36 | 10.2 | 0.0709 | 0.0709 | 2.19 | 2.07 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | paralog_exclusion | 17685089 | 17685 | 3.82 | 2.92 | 4.77 | 7.28e-36 | 1.34 | 10.9 | 0.0689 | 0.0649 | 2.25 | 2.09 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_i | 7038496 | 7.04e+03 | 1.61 | 1.01 | 2.29 | 0.0013 | 0.78 | 3.34 | 0.111 | 0.135 | 2 | 1.84 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_iii | 10033495 | 10033 | 0.227 | 2.28 | 0.316 | 0.988 | 0.0179 | 2.89 | 0.516 | 0.744 | 0.955 | 0.482 |
| control | M1 |  | hct116 | TEST | NC2_M1 | 8697840 | 8.7e+03 | 1.46 | 2.43 | 0.678 | 0.233 | 0.243 | 8.75 | 0.0679 | 0.0589 | 2.85 | 2.87 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | NC2 | 8697840 | 8.7e+03 | 0.121 | 1.21 | 0.241 | 1 | 0.0172 | 0.853 | 0.754 | 0.647 | 0.494 | 0.949 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | no_split_sibling | 7074275 | 7.07e+03 | 2.61 | 1.63 | 3.16 | 2.02e-06 | 0.411 | 16.6 | 0.151 | 0.155 | 2.42 | 2.39 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | TEST |  | 17692753 | 17693 | 3.76 | 3.81 | 4.65 | 1.14e-34 | 2.04 | 6.94 | 0.0839 | 0.0819 | 2.07 | 1.99 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | TEST |  | 17692753 | 17693 | 17.3 | 20.3 | 12.9 | 0 | 10.6 | 28.2 | 0.022 | 0.027 | 1.71 | 1.53 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.492 | 1.65 | 0.63 | 0.982 | 0.0173 | 14.1 | 0.878 | 0.793 | 0.238 | 0.436 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.492 | 0.41 | 0.531 | 0.982 | 0.126 | 1.92 | 0.782 | 0.845 | 0.577 | 0.576 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | TEST |  | 8373505 | 8.37e+03 | 5.06 | 15.1 | 6.35 | 1.58e-23 | 0.881 | 29.1 | 0.905 | 0.706 | 0.352 | 0.58 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.082 | 0.41 | 0.253 | 1 | 0.0187 | 0.359 | 0.982 | 0.669 | 0.0308 | 0.481 |
| baseline | B0 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.739 | 0.41 | 0.853 | 0.856 | 0.414 | 1.32 | 0.754 | 0.561 | 0.789 | 0.95 |
| baseline | B1 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.986 | 0.41 | 1.24 | 0.558 | 0.253 | 3.84 | 0.615 | 0.641 | 0.825 | 0.874 |
| baseline | B2 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 1.4 | 0.41 | 2.41 | 0.11 | 0.255 | 7.66 | 0.267 | 0.517 | 1.37 | 0.956 |
| baseline | B3 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.739 | 0.41 | 7.78 | 0.856 | 7.26e-07 | 752532 | 0.285 | 0.05 | 4.35 | 2.78 |
| baseline | B4 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.986 | 0.41 | 1.14 | 0.558 | 0.441 | 2.2 | 0.965 | 0.97 | 0.834 | 0.827 |
| baseline | B4e |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 14.5 | 17.7 | 17.9 | 1.31e-130 | 10.2 | 20.7 | 0.485 | 0.485 | 1.01 | 1.01 |
| baseline | B5 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 35.9 | 4.96 | 13.4 | 0 | 5.28 | 244 | 0.145 | 0.152 | 1.92 | 1.56 |
| baseline | B6 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.41 | 1.65 | 0.804 | 0.993 | 0.0592 | 2.84 | 0.784 | 0.688 | 0.681 | 0.795 |
| control | NC1 | NC1: residual-permutation null | hek293t | TEST |  | 8373505 | 8.37e+03 | 1.73 | 0.822 | 1.29 | 0.0132 | 0.2 | 14.9 | 0.107 | 0.163 | 1.69 | 1.43 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | TEST |  | 8373505 | 8.37e+03 | 12.8 | 26.5 | 7.29 | 6.9e-108 | 5.87 | 27.8 | 0.34 | 0.13 | 1.13 | 1.34 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 1.48 | 0.822 | 1.34 | 0.0694 | 0.222 | 9.9 | 0.573 | 0.459 | 0.82 | 1.22 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.492 | 1.65 | 0.622 | 0.982 | 0.0173 | 14.1 | 0.878 | 0.788 | 0.238 | 0.413 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.575 | 0.822 | 0.713 | 0.958 | 0.0134 | 24.7 | 0.86 | 0.674 | 0.28 | 0.663 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.492 | 1.65 | 0.605 | 0.982 | 0.0158 | 15.4 | 0.844 | 0.768 | 0.288 | 0.442 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.492 | 1.65 | 0.63 | 0.982 | 0.0173 | 14 | 0.88 | 0.793 | 0.238 | 0.415 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | paralog_exclusion | 8369889 | 8.37e+03 | 0.493 | 1.65 | 0.63 | 0.982 | 0.0173 | 14.1 | 0.879 | 0.796 | 0.238 | 0.435 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_i | 755076 | 755 | 0.495 | 1.22 | 0.172 | 0.911 | 0.0164 | 14.9 | 0.495 | 0.409 | 0.975 | 1.9 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_iii | 7524519 | 7.52e+03 | 0.0615 | 0.616 | 0.0366 | 1 | 0.00874 | 0.433 | 0.763 | 0.56 | 0.483 | 0.948 |
| control | M1 |  | hek293t | TEST | NC2_M1 | 6331578 | 6.33e+03 | 0.324 | 1.62 | 1.2 | 0.985 | 0.0445 | 2.35 | 0.799 | 0.641 | 0.477 | 0.945 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | NC2 | 6331578 | 6.33e+03 | 0.0809 | 0.809 | 0.0963 | 1 | 8.74e-05 | 74.8 | 0.655 | 0.507 | 0.95 | 0.95 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | no_split_sibling | 3371874 | 3.37e+03 | 0.114 | 1.14 | 0.407 | 1 | 0.0395 | 0.327 | 0.92 | 0.879 | 0.122 | 0.247 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 0.739 | 0.41 | 0.597 | 0.856 | 0.0491 | 11.1 | 0.853 | 0.743 | 0.304 | 0.502 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | TEST |  | 8373505 | 8.37e+03 | 4.73 | 16.8 | 6.38 | 9.54e-21 | 0.684 | 32.6 | 0.917 | 0.711 | 0.339 | 0.582 |
| control | PC2 |  | hct116 | TEST |  | 5822281 | 5.82e+03 | 181 | 1.02e+03 | 20.2 | 0 | 75.2 | 435 | 0.000999 | 0.000999 | 59.1 | 7.08 |

## cmh

| role | run | label | line | split | subset | CMH_OR(reg x S2 expr) | CMH_OR(own class x S2 expr) | CMH_OR(reg x HCT116 expr) [D6 diagnostic] | mean_AUROC | top10_share | OR_hub_excluded |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST |  | 0.942 | 1.8 | 0.993 | 0.5 | 0.261 | 2.69 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | TEST |  | 0.909 | 1.24 | 0.973 | 0.498 | 0.161 | 1.39 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | TEST |  | 3.4 | 3.44 | 2.82 | 0.646 | 0.123 | 14.8 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | TEST |  | 1.2 | 2.01 | 1.2 | 0.504 | 0.291 | 2.19 |
| baseline | B0 |  | hct116 | TEST |  | 1.04 | 1.18 | 1.03 | 0.5 | 0.0109 | 0.993 |
| baseline | B1 |  | hct116 | TEST |  | 3.13 | 1.4 | 2.26 | 0.514 | 0.133 | 1.39 |
| baseline | B2 |  | hct116 | TEST |  | 0.958 | 3.08 | 0.808 | 0.482 | 0.227 | 3.99 |
| baseline | B3 |  | hct116 | TEST |  |  | 0.0133 |  | 0.5 | 1 | 0.0992 |
| baseline | B4 |  | hct116 | TEST |  | 7.15 | 4.6 | 2.9 | 0.577 | 0.00622 | 2.79 |
| baseline | B4e |  | hct116 | TEST |  | 4.12 | 2.73 | 3.62 | 0.699 | 0.00622 | 8.03 |
| baseline | B5 |  | hct116 | TEST |  | 3.96 | 4.22 | 2.92 | 0.699 | 0.387 | 11.8 |
| baseline | B6 |  | hct116 | TEST |  | 0.687 | 0.563 | 0.747 | 0.506 | 0.017 | 0.695 |
| control | NC1 | NC1: residual-permutation null | hct116 | TEST |  | 1.11 | 1.32 | 1.16 | 0.484 | 0.0994 | 0.894 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | TEST |  | 2.49 | 2.13 | 2.03 | 0.633 | 0.0711 | 9.97 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | TEST |  | 0.731 | 1.6 | 0.799 | 0.504 | 0.335 | 1.89 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | TEST |  | 0.928 | 1.77 | 0.98 | 0.5 | 0.262 | 2.69 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | TEST |  | 0.82 | 1.29 | 0.844 | 0.498 | 0.286 | 2.29 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | TEST |  | 0.957 | 1.64 | 1.01 | 0.501 | 0.267 | 2.59 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | TEST |  | 0.918 | 1.74 | 0.97 | 0.5 | 0.262 | 2.59 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | paralog_exclusion | 0.932 | 1.79 | 0.983 | 0.5 | 0.261 | 2.69 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_i | 0.841 | 1.52 | 0.904 | 0.5 | 0.305 | 0.876 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_iii | 0.495 | 0.216 | 0.497 | 0.499 | 0.997 | 0.45 |
| control | M1 |  | hct116 | TEST | NC2_M1 | 2.2 | 1.45 | 2.52 | 0.496 | 0.102 | 2.94 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | NC2 | 0 | 0 | 0 | 0.5 | 1 | 0.487 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | no_split_sibling | 0.64 | 1.44 | 0.677 | 0.491 | 0.459 | 1.06 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | TEST |  | 0.977 | 1.82 | 0.991 | 0.497 | 0.263 | 1.99 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | TEST |  | 2.98 | 3.25 | 2.44 | 0.646 | 0.125 | 12.6 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST |  | 0.675 | 0.139 | 0.683 | 0.5 | 0.796 | 0.833 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | TEST |  | 0.961 | 0.441 | 1.05 | 0.507 | 0.255 | 0.625 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | TEST |  | 3.11 | 0.515 | 3.36 | 0.517 | 0.537 | 3.98 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | TEST |  | 0.24 | 0.0279 | 0.25 | 0.498 | 0.874 | 0.208 |
| baseline | B0 |  | hek293t | TEST |  | 0.797 | 0.734 | 0.796 | 0.498 | 0.0233 | 0.833 |
| baseline | B1 |  | hek293t | TEST |  | 1.62 | 0.885 | 1.59 | 0.532 | 0.239 | 1.04 |
| baseline | B2 |  | hek293t | TEST |  | 0.926 | 0.88 | 0.908 | 0.497 | 0.349 | 1.25 |
| baseline | B3 |  | hek293t | TEST |  | 1.12 | 0.198 | 1.17 | 0.5 | 1 | 1.88 |
| baseline | B4 |  | hek293t | TEST |  | 0.591 | 0.601 | 0.44 | 0.453 | 0.0131 | 1.04 |
| baseline | B4e |  | hek293t | TEST |  | 5.62 | 3.76 | 5.89 | 0.691 | 0.0143 | 11.7 |
| baseline | B5 |  | hek293t | TEST |  | 2.54 | 4.91 | 2.86 | 0.691 | 0.984 | 15.3 |
| baseline | B6 |  | hek293t | TEST |  | 0.482 | 0.403 | 0.5 | 0.506 | 0.0406 | 0.104 |
| control | NC1 | NC1: residual-permutation null | hek293t | TEST |  | 1.01 | 1.91 | 1.04 | 0.497 | 0.167 | 0.416 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | TEST |  | 1.32 | 1.65 | 1.51 | 0.514 | 0.253 | 5.03 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | TEST |  | 1.13 | 0.488 | 1.17 | 0.499 | 0.822 | 2.51 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | TEST |  | 0.672 | 0.139 | 0.68 | 0.5 | 0.798 | 0.833 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | TEST |  | 0.799 | 0.15 | 0.827 | 0.499 | 0.791 | 1.04 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | TEST |  | 0.732 | 0.158 | 0.75 | 0.5 | 0.807 | 0.625 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | TEST |  | 0.673 | 0.139 | 0.68 | 0.5 | 0.798 | 0.833 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | paralog_exclusion | 0.675 | 0.139 | 0.682 | 0.5 | 0.796 | 0.833 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_i | 3.49 | 0.452 | 3.47 | 0.498 | 0.915 | 0.399 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_iii | 0 | 0 | 0 | 0.5 | 1 | 0.168 |
| control | M1 |  | hek293t | TEST | NC2_M1 | 0.547 | 0.334 | 0.558 | 0.505 | 0.132 | 0.232 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | NC2 | 0 | 0 | 0 | 0.5 | 1 | 0.232 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | no_split_sibling | 0 | 0 | 0 | 0.501 | 0.973 | 0.365 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | TEST |  | 1.08 | 0.213 | 1.05 | 0.498 | 0.801 | 1.25 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | TEST |  | 2.78 | 0.479 | 3.05 | 0.516 | 0.529 | 4.4 |
| control | PC2 |  | hct116 | TEST |  | 7.55 | 65.1 | 7.35 | 0.556 | 0.666 | 169 |

## sign

| role | run | label | line | split | subset | sign_concordant_OR | n_overlap | concordance | c0 | sign_OR | sign_perm_p | reg_weighted_concordance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST |  | 5.84 | 131 | 0.718 | 0.541 | 7.02 | 0.012 | 0.654 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | TEST |  | 1.89 | 50 | 0.72 | 0.493 | 12.8 | 0.00599 | 0.902 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | TEST |  | 21.1 | 601 | 0.549 | 0.507 | 1.45 | 0.216 | 0.491 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | TEST |  | 7.18 | 153 | 0.634 | 0.607 | 1.41 | 0.905 | 0.516 |
| baseline | B0 |  | hct116 | TEST |  | 1.25 | 40 | 0.7 | 0.7 | 2.28 | 1 | 0.714 |
| baseline | B1 |  | hct116 | TEST |  | 0.989 | 39 | 0.462 | 0.547 | 0.409 | 1 | 0.545 |
| baseline | B2 |  | hct116 | TEST |  | 7.83 | 204 | 0.843 | 0.845 | 1.83 | 1 | 0.82 |
| baseline | B3 |  | hct116 | TEST |  | 0.0222 | 1 | 0 | 0 | 0.333 | 1 |  |
| baseline | B4 |  | hct116 | TEST |  | 3.68 | 194 | 0.423 | 0.423 | 0.733 | 1 | 0.256 |
| baseline | B4e |  | hct116 | TEST |  | 5.9 | 249 | 0.526 | 0.526 | 1.11 | 1 | 0.548 |
| baseline | B5 |  | hct116 | TEST |  | 26.2 | 736 | 0.757 | 0.757 | 3.11 | 1 | 0.702 |
| baseline | B6 |  | hct116 | TEST |  | 0.585 | 21 | 0.524 | 0.524 | 1 | 0.465 |  |
| control | NC1 | NC1: residual-permutation null | hct116 | TEST |  | 0.828 | 44 | 0.364 | 0.417 | 0.548 | 0.517 | 0.0952 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | TEST |  | 13 | 354 | 0.585 | 0.514 | 1.82 | 0.000999 | 0.616 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | TEST |  | 4.55 | 116 | 0.586 | 0.467 | 3.29 | 0.133 | 0.524 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | TEST |  | 5.73 | 129 | 0.713 | 0.528 | 7.11 | 0.004 | 0.656 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | TEST |  | 4.03 | 99 | 0.657 | 0.497 | 7.5 | 0.012 | 0.605 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | TEST |  | 5.09 | 119 | 0.689 | 0.523 | 5.65 | 0.00699 | 0.655 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | TEST |  | 5.6 | 127 | 0.709 | 0.534 | 6.6 | 0.005 | 0.651 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | paralog_exclusion | 5.79 | 130 | 0.715 | 0.539 | 6.91 | 0.004 | 0.654 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_i | 2.3 | 48 | 0.771 | 0.521 | 13 | 0.00799 | 0.889 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_iii | 13 | 1 | 0 | 0 | 0.333 | 1 |  |
| control | M1 |  | hct116 | TEST | NC2_M1 | 1.6 | 6 | 0.5 | 0.5 | 1 | 1 |  |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | NC2 | 47.6 | 0 |  |  |  |  |  |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | no_split_sibling | 3.59 | 32 | 0.5 | 0.477 | 1.25 | 0.481 | 0.565 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | TEST |  | 3.71 | 128 | 0.414 | 0.469 | 0.575 | 0.996 | 0.336 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | TEST |  | 20.6 | 567 | 0.504 | 0.488 | 1.15 | 0.526 | 0.47 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST |  | 1.03 | 6 | 0.333 | 0.556 | 0.2 | 1 | 0.2 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | TEST |  | 0.156 | 6 | 0.167 | 0.5 | 0.0857 | 1 |  |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | TEST |  | 12 | 61 | 0.377 | 0.404 | 0.306 | 0.383 | 0.645 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | TEST |  | 0.306 | 1 | 0 | 0 | 0.333 | 1 |  |
| baseline | B0 |  | hek293t | TEST |  | 0.325 | 9 | 0.222 | 0.222 | 0.333 | 1 |  |
| baseline | B1 |  | hek293t | TEST |  | 1.29 | 12 | 0.667 | 0.5 | 7.22 | 1 |  |
| baseline | B2 |  | hek293t | TEST |  | 1.62 | 17 | 0.588 | 0.588 | 1.4 | 1 | 0.8 |
| baseline | B3 |  | hek293t | TEST |  | 0.0832 | 9 | 0 | 0 | 0.0526 | 1 | 0 |
| baseline | B4 |  | hek293t | TEST |  | 0.649 | 12 | 0.333 | 0.333 | 0.529 | 1 |  |
| baseline | B4e |  | hek293t | TEST |  | 11.9 | 171 | 0.415 | 0.415 | 0.711 | 1 | 0.462 |
| baseline | B5 |  | hek293t | TEST |  | 30.9 | 403 | 0.449 | 0.449 | 0.816 | 1 | 0.324 |
| baseline | B6 |  | hek293t | TEST |  | 0.152 | 5 | 0.2 | 0.2 | 0.333 | 1 |  |
| control | NC1 | NC1: residual-permutation null | hek293t | TEST |  | 0.326 | 21 | 0.0952 | 0.308 | 0.0208 | 1 | 0 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | TEST |  | 56.5 | 151 | 0.57 | 0.525 | 1.62 | 0.987 | 0.532 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | TEST |  | 2.67 | 18 | 0.278 | 0.463 | 0.0748 | 1 | 0.278 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | TEST |  | 1.04 | 6 | 0.333 | 0.556 | 0.2 | 1 | 0.2 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | TEST |  | 1.02 | 7 | 0.286 | 0.531 | 0.143 | 1 | 0.167 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | TEST |  | 1.04 | 6 | 0.333 | 0.5 | 0.238 | 1 | 0.2 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | TEST |  | 1.04 | 6 | 0.333 | 0.556 | 0.2 | 1 | 0.2 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | paralog_exclusion | 1.03 | 6 | 0.333 | 0.556 | 0.2 | 1 | 0.2 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_i | 0.272 | 2 | 0 | 0.5 | 0.111 | 1 |  |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_iii | 4.87 | 0 |  |  |  |  |  |
| control | M1 |  | hek293t | TEST | NC2_M1 | 0.622 | 2 | 1 | 0.5 | 9 | 1 |  |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | NC2 | 4.97 | 0 |  |  |  |  |  |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | no_split_sibling | 0.81 | 0 |  |  |  |  |  |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | TEST |  | 1.1 | 9 | 0.222 | 0.346 | 0.167 | 1 | 0.143 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | TEST |  | 14 | 57 | 0.439 | 0.492 | 0.123 | 1 | 0.621 |
| control | PC2 |  | hct116 | TEST |  | 9.23e+03 | 301 | 0.99 | 0.703 | 7.51e+03 | 0.000999 | 0.997 |

## callset

| role | run | label | line | split | subset | n_calls | callset_OR | callset_fisher_p |
|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST |  | 1187033 | 1.15 | 3.28e-11 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | TEST |  | 1187033 | 1.15 | 3.28e-11 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | TEST |  | 1187033 | 1.15 | 3.28e-11 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | TEST |  | 0 | 514 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | TEST |  | 367013 | 1.11 | 0.00315 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | TEST |  | 1186476 | 1.14 | 4.16e-10 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | TEST |  | 3226477 | 1.13 | 1.51e-20 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | TEST |  | 1099644 | 1.12 | 5.01e-08 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | TEST |  | 1181215 | 1.14 | 1.69e-10 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | paralog_exclusion | 1186551 | 1.15 | 3.34e-11 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_i | 524353 | 1.08 | 0.00025 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_iii | 621137 | 0.764 | 1 |
| control | M1 |  | hct116 | TEST | NC2_M1 | 529899 | 0.803 | 0.999 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | NC2 | 529899 | 0.803 | 0.999 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | no_split_sibling | 646760 | 0.884 | 1 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST |  | 512822 | 1.16 | 2.16e-05 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | TEST |  | 512822 | 1.16 | 2.16e-05 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | TEST |  | 512822 | 1.16 | 2.16e-05 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | TEST |  | 0 | 687 | 1 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | TEST |  | 150101 | 1.2 | 0.00237 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | TEST |  | 512743 | 1.16 | 2.11e-05 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | TEST |  | 1417212 | 1.14 | 1.71e-08 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | TEST |  | 474843 | 1.15 | 0.000137 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | TEST |  | 510544 | 1.16 | 1.69e-05 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | paralog_exclusion | 512629 | 1.16 | 2.15e-05 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_i | 56025 | 0.786 | 1 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_iii | 451020 | 1.29 | 3.53e-09 |
| control | M1 |  | hek293t | TEST | NC2_M1 | 379315 | 1.18 | 0.000487 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | NC2 | 379315 | 1.18 | 0.000487 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | no_split_sibling | 269793 | 1.12 | 0.0153 |
| control | PC2 |  | hct116 | TEST |  | 9.89e+03 | 119 | 0 |

## C1

| role | run | label | line | split | subset | jk_CI95_lo | jk_CI95_hi | dlogOR | OR_run | OR_C1 | dOR | c1_mode |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| control | g1_t100_w0__a3_r095_s2_w0 | C1 pair-specific gain | hct116 | TEST |  | -0.747 | 0.794 | 0.0234 | 3.85 | 3.76 | 0.0891 | score (approximation) |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | C1 pair-specific gain | hct116 | TEST |  | -0.0635 | 0.186 | 0.0612 | 18.4 | 17.3 | 1.09 | pair_score x target_weight |
| control | g1_t100_w0__a3_r095_s2_w0 | C1 pair-specific gain | hek293t | TEST |  | -1.09 | 0.278 | -0.406 | 0.492 | 0.739 | -0.247 | score (approximation) |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | C1 pair-specific gain | hek293t | TEST |  | -0.143 | 0.28 | 0.0686 | 5.06 | 4.73 | 0.336 | pair_score x target_weight |

## paired

| role | run | label | line | split | subset | jk_CI95_lo | jk_CI95_hi | dlogOR | boot_lb90 | jk_lb90 | dlogCMH | boot_cmh_ub90 | boot_cmh_n_undefined | adopt_style_no_folds | adopt_rule_8_2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs a3_r095_s2_w0 | hct116 | TEST |  | 0.297 | 1.64 | 0.97 | 0.449 | 0.418 | 0.035 | 0.423 | 0 | True | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | hct116 | TEST |  | -2.42 | -0.705 | -1.56 | -2.21 | -2.27 | -1.28 | -0.697 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs focal__g1_t100_w0__a3_r095_s2_w0 | hct116 | TEST |  | -0.688 | 0.374 | -0.157 | -0.454 | -0.592 | -0.24 | 0.492 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B0 | hct116 | TEST |  | 0.13 | 2.26 | 1.19 | 0.237 | 0.322 | -0.101 | 0.689 | 0 | True | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B1 | hct116 | TEST |  | -0.478 | 2.92 | 1.22 | -0.00334 | -0.172 | -1.2 | -0.116 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B2 | hct116 | TEST |  | -1.69 | 0.787 | -0.449 | -1.6 | -1.46 | -0.0168 | 0.517 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B3 | hct116 | TEST |  | -14.7 | 24.5 | 4.89 | -5.53 | -11.2 |  |  | 592 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B4 | hct116 | TEST |  | -1.7 | 0.906 | -0.398 | -1.49 | -1.47 | -2.03 | -1.24 | 2 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B4e | hct116 | TEST |  | -1.93 | 0.623 | -0.652 | -1.66 | -1.7 | -1.48 | -0.853 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B5 | hct116 | TEST |  | -4.13 | 0.57 | -1.78 | -2.97 | -3.7 | -1.44 | -0.813 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B6 | hct116 | TEST |  | 0.5 | 3.18 | 1.84 | 0.779 | 0.741 | 0.316 | 1.16 | 0 | True | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs NC1 | hct116 | TEST |  | -0.656 | 2.85 | 1.1 | 0.0488 | -0.34 | -0.16 | 0.763 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs C2_g2 | hct116 | TEST |  | -1.99 | -0.0374 | -1.01 | -1.8 | -1.81 | -0.971 | -0.28 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__S1 | hct116 | TEST |  | -0.291 | 0.536 | 0.123 | -0.317 | -0.216 | 0.253 | 0.634 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__noHCT116 | hct116 | TEST |  | -0.0474 | 0.0785 | 0.0156 | -0.0361 | -0.036 | 0.0146 | 0.073 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__hidden_factors | hct116 | TEST |  | -0.0732 | 0.639 | 0.283 | -0.0295 | -0.00905 | 0.139 | 0.422 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion | hct116 | TEST |  | -0.0553 | 0.25 | 0.0971 | -0.0579 | -0.0278 | -0.0154 | 0.102 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop | hct116 | TEST |  | -0.0297 | 0.0924 | 0.0314 | 0.00265 | -0.0187 | 0.0258 | 0.0633 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 vs B5 | hct116 | TEST |  | -2.29 | 1.86 | -0.215 | -1.19 | -1.91 | -0.151 | 0.122 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 vs C2_g2 | hct116 | TEST |  | 0.0717 | 1.03 | 0.551 | 0.111 | 0.158 | 0.314 | 0.683 | 0 | True | False |
| comparison | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 vs g1_t100_w0__a3_r095_s2_w0 | hct116 | TEST |  | 0.705 | 2.42 | 1.56 | 1.06 | 0.86 | 1.28 | 1.85 | 0 | True | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs a3_r095_s2_w0 | hek293t | TEST |  | -3.22 | 3.22 | 0 | -2.66 | -2.64 | -0.353 |  | 45 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | hek293t | TEST |  | -4.2 | -0.464 | -2.33 | -4.36 | -3.86 | -1.53 | -0.646 | 6 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs focal__g1_t100_w0__a3_r095_s2_w0 | hek293t | TEST |  | -2.13 | 5.72 | 1.79 | -0.998 | -1.42 | 1.03 |  | 150 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B0 | hek293t | TEST |  | -3.91 | 3.1 | -0.406 | -3.02 | -3.28 | -0.166 | 1.47 | 9 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B1 | hek293t | TEST |  | -4.75 | 3.36 | -0.694 | -3.6 | -4.02 | -0.877 | 1.46 | 2 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B2 | hek293t | TEST |  | -4.94 | 2.85 | -1.04 | -3.85 | -4.24 | -0.316 | 0.76 | 19 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B3 | hek293t | TEST |  | -14.2 | 13.3 | -0.406 | -9.08 | -11.7 | -0.504 | 0.535 | 330 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B4 | hek293t | TEST |  | -4.06 | 2.67 | -0.694 | -3.34 | -3.45 | 0.134 | 1.86 | 10 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B4e | hek293t | TEST |  | -6.75 | -0.0123 | -3.38 | -5.97 | -6.15 | -2.12 | -0.867 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B5 | hek293t | TEST |  | -8.55 | -0.0292 | -4.29 | -7.38 | -7.78 | -1.32 | -0.0824 | 7 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs B6 | hek293t | TEST |  | -4.03 | 4.39 | 0.183 | -2.72 | -3.27 | 0.337 |  | 35 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs NC1 | hek293t | TEST |  | -5.23 | 2.72 | -1.26 | -3.78 | -4.51 | -0.402 | 0.995 | 13 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs C2_g2 | hek293t | TEST |  | -6.88 | 0.368 | -3.26 | -5.94 | -6.22 | -0.674 | 0.405 | 1 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__S1 | hek293t | TEST |  | -3.72 | 1.52 | -1.1 | -3.35 | -3.25 | -0.512 |  | 53 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__noHCT116 | hek293t | TEST |  | -0.00245 | 0.00245 | 0 | -0.00388 | -0.00201 | 0.00509 | 0.0122 | 197 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__hidden_factors | hek293t | TEST |  | -0.61 | 0.301 | -0.154 | -0.511 | -0.528 | -0.169 | 0.168 | 196 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion | hek293t | TEST |  | -0.325 | 0.325 | 0 | -0.496 | -0.266 | -0.0809 | 0.42 | 195 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0 | g1_t100_w0__a3_r095_s2_w0 vs g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop | hek293t | TEST |  | -0.00581 | 0.00581 | 0 | -0.00732 | -0.00476 | 0.00295 | 0.007 | 197 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 vs B5 | hek293t | TEST |  | -4.48 | 0.56 | -1.96 | -4.1 | -4.02 | 0.206 | 1.15 | 7 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 vs C2_g2 | hek293t | TEST |  | -2.95 | 1.1 | -0.925 | -3.13 | -2.59 | 0.855 | 1.59 | 0 | False | False |
| comparison | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 vs g1_t100_w0__a3_r095_s2_w0 | hek293t | TEST |  | 0.464 | 4.2 | 2.33 | 1.16 | 0.8 | 1.53 |  | 6 | True | False |

## runtime

| role | run | label | line | split | subset | eval_wall_s | score_wall_s | threads |
|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST |  | 2.57 | 67.2 | {'cpu_count': 10, 'blas': []} |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | TEST |  | 2.54 | 66.1 | {'cpu_count': 10, 'blas': []} |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | TEST |  | 3.2 | 82.7 | {'cpu_count': 10, 'blas': []} |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | TEST |  | 2.44 | 51.5 | {'cpu_count': 10, 'blas': []} |
| baseline | B0 |  | hct116 | TEST |  | 2.24 | 2.73 | {'cpu_count': 10, 'blas': []} |
| baseline | B1 |  | hct116 | TEST |  | 2.17 | 4.12 | {'cpu_count': 10, 'blas': []} |
| baseline | B2 |  | hct116 | TEST |  | 2.19 | 4.8 | {'cpu_count': 10, 'blas': []} |
| baseline | B3 |  | hct116 | TEST |  | 1.73 | 2.39 | {'cpu_count': 10, 'blas': []} |
| baseline | B4 |  | hct116 | TEST |  | 2.05 | 2.68 | {'cpu_count': 10, 'blas': []} |
| baseline | B4e |  | hct116 | TEST |  | 2.2 | 2.86 | {'cpu_count': 10, 'blas': []} |
| baseline | B5 |  | hct116 | TEST |  | 1.98 | 2.41 | {'cpu_count': 10, 'blas': []} |
| baseline | B6 |  | hct116 | TEST |  | 1.92 | 12.2 | {'cpu_count': 10, 'blas': []} |
| control | NC1 | NC1: residual-permutation null | hct116 | TEST |  | 2.26 | 27.5 | {'cpu_count': 10, 'blas': []} |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | TEST |  | 1.92 | 42.8 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | TEST |  | 2.52 | 53.9 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | TEST |  | 2.58 | 66.5 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | TEST |  | 2.62 | 90.7 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | TEST |  | 2.57 | 66.6 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | TEST |  | 2.52 | 66.1 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | paralog_exclusion | 2.65 | 67.2 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_i | 1.82 | 67.2 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_iii | 1.8 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | M1 |  | hct116 | TEST | NC2_M1 | 1.86 | 51.6 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | NC2 | 1.73 | 67.2 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | no_split_sibling | 1.78 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | TEST |  | 0.963 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | TEST |  | 1.07 | 82.7 | {'cpu_count': 10, 'blas': []} |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST |  | 1.78 | 67.2 | {'cpu_count': 10, 'blas': []} |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | TEST |  | 2.14 | 66.1 | {'cpu_count': 10, 'blas': []} |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | TEST |  | 2.24 | 82.7 | {'cpu_count': 10, 'blas': []} |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | TEST |  | 2.01 | 51.5 | {'cpu_count': 10, 'blas': []} |
| baseline | B0 |  | hek293t | TEST |  | 1.62 | 2.73 | {'cpu_count': 10, 'blas': []} |
| baseline | B1 |  | hek293t | TEST |  | 1.53 | 4.12 | {'cpu_count': 10, 'blas': []} |
| baseline | B2 |  | hek293t | TEST |  | 1.69 | 4.8 | {'cpu_count': 10, 'blas': []} |
| baseline | B3 |  | hek293t | TEST |  | 1.42 | 2.39 | {'cpu_count': 10, 'blas': []} |
| baseline | B4 |  | hek293t | TEST |  | 1.57 | 2.68 | {'cpu_count': 10, 'blas': []} |
| baseline | B4e |  | hek293t | TEST |  | 1.66 | 2.86 | {'cpu_count': 10, 'blas': []} |
| baseline | B5 |  | hek293t | TEST |  | 1.57 | 2.41 | {'cpu_count': 10, 'blas': []} |
| baseline | B6 |  | hek293t | TEST |  | 1.62 | 12.2 | {'cpu_count': 10, 'blas': []} |
| control | NC1 | NC1: residual-permutation null | hek293t | TEST |  | 1.7 | 27.5 | {'cpu_count': 10, 'blas': []} |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | TEST |  | 1.45 | 42.8 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | TEST |  | 2.01 | 53.9 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | TEST |  | 2.02 | 66.5 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | TEST |  | 2.07 | 90.7 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | TEST |  | 2.01 | 66.6 | {'cpu_count': 10, 'blas': []} |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | TEST |  | 2.03 | 66.1 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | paralog_exclusion | 2.08 | 67.2 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_i | 1.47 | 67.2 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_iii | 1.72 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | M1 |  | hek293t | TEST | NC2_M1 | 1.68 | 51.6 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | NC2 | 1.66 | 67.2 | {'cpu_count': 10, 'blas': []} |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | no_split_sibling | 1.59 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | TEST |  | 0.479 | 67.2 | {'cpu_count': 10, 'blas': []} |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | TEST |  | 0.472 | 82.7 | {'cpu_count': 10, 'blas': []} |
| control | PC2 |  | hct116 | TEST |  | 1.51 |  |  |

## class:RBP

| role | run | label | line | split | subset | OR@0.1% | fisher_p(indep.) | n_regulators | selected | positives | a |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST |  | 3.56 | 1.21e-08 | 17 | 482 | 3.06e+03 | 30 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | TEST |  | 4.48 | 9.05e-07 | 17 | 220 | 3.06e+03 | 17 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | TEST |  | 22 | 6.09e-66 | 17 | 261 | 3.06e+03 | 75 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | TEST |  | 5.25 | 3.27e-10 | 17 | 269 | 3.06e+03 | 24 |
| baseline | B0 |  | hct116 | TEST |  | 1.46 | 0.299 | 17 | 150 | 3.06e+03 | 4 |
| baseline | B1 |  | hct116 | TEST |  | 0.131 | 1 | 17 | 202 | 3.06e+03 | 0 |
| baseline | B2 |  | hct116 | TEST |  | 2.69 | 6.08e-08 | 17 | 860 | 3.06e+03 | 41 |
| baseline | B3 |  | hct116 | TEST |  | 53.3 | 1 | 17 | 0 | 3.06e+03 | 0 |
| baseline | B4 |  | hct116 | TEST |  | 2.38 | 0.0331 | 17 | 164 | 3.06e+03 | 7 |
| baseline | B4e |  | hct116 | TEST |  | 4.93 | 2.92e-06 | 17 | 166 | 3.06e+03 | 14 |
| baseline | B5 |  | hct116 | TEST |  | 0.143 | 1 | 17 | 186 | 3.06e+03 | 0 |
| baseline | B6 |  | hct116 | TEST |  | 0.745 | 0.749 | 17 | 145 | 3.06e+03 | 2 |
| control | NC1 | NC1: residual-permutation null | hct116 | TEST |  | 0.832 | 0.701 | 17 | 65 | 3.06e+03 | 1 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | TEST |  | 6.65 | 6.27e-10 | 17 | 172 | 3.06e+03 | 19 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | TEST |  | 3.06 | 4.65e-08 | 17 | 631 | 3.06e+03 | 34 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | TEST |  | 3.65 | 4.14e-09 | 17 | 487 | 3.06e+03 | 31 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | TEST |  | 3.35 | 2.18e-09 | 17 | 614 | 3.06e+03 | 36 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | TEST |  | 3.77 | 3.62e-09 | 17 | 457 | 3.06e+03 | 30 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | TEST |  | 3.52 | 1.53e-08 | 17 | 487 | 3.06e+03 | 30 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | paralog_exclusion | 3.57 | 1.17e-08 | 17 | 482 | 3.05e+03 | 30 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_i | 2.76 | 0.000232 | 10 | 221 | 3.05e+03 | 18 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_iii | 13061 | 1 | 6 | 0 | 4 | 0 |
| control | M1 |  | hct116 | TEST | NC2_M1 | 221 | 1 | 6 | 29 | 4 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | NC2 | 13061 | 1 | 6 | 0 | 4 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | no_split_sibling | 1.45 | 0.0604 | 3 | 334 | 1.43e+03 | 23 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | TEST |  | 5.67 | 1.13e-19 | 17 | 494 | 3.06e+03 | 47 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | TEST |  | 18.3 | 1.77e-52 | 17 | 254 | 3.06e+03 | 64 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | TEST |  | 10.3 | 1 | 4 | 7 | 267 | 0 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| baseline | B0 |  | hek293t | TEST |  | 4.87 | 0.191 | 4 | 33 | 267 | 1 |
| baseline | B1 |  | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| baseline | B2 |  | hek293t | TEST |  | 9.12 | 1 | 4 | 8 | 267 | 0 |
| baseline | B3 |  | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| baseline | B4 |  | hek293t | TEST |  | 1.91 | 1 | 4 | 40 | 267 | 0 |
| baseline | B4e |  | hek293t | TEST |  | 12.4 | 0.0023 | 4 | 41 | 267 | 3 |
| baseline | B5 |  | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| baseline | B6 |  | hek293t | TEST |  | 1.78 | 1 | 4 | 43 | 267 | 0 |
| control | NC1 | NC1: residual-permutation null | hek293t | TEST |  | 8.16 | 1 | 4 | 9 | 267 | 0 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | TEST |  | 87.4 | 5.41e-56 | 4 | 120 | 267 | 39 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | TEST |  | 0.874 | 1 | 4 | 88 | 267 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | paralog_exclusion | 155 | 1 | 4 | 0 | 267 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_i | 76.6 | 1 | 2 | 0 | 267 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_iii | 41969 | 1 | 2 | 0 | 0 | 0 |
| control | M1 |  | hek293t | TEST | NC2_M1 | 791 | 1 | 2 | 26 | 0 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | NC2 | 41969 | 1 | 2 | 0 | 0 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | TEST |  | 155 | 1 | 4 | 0 | 267 | 0 |
| control | PC2 |  | hct116 | TEST |  | 170 | 8.72e-21 | 3 | 198 | 23 | 12 |

## class:TF

| role | run | label | line | split | subset | OR@0.1% | fisher_p(indep.) | n_regulators | selected | positives | a |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST |  | 6.67 | 1 | 137 | 82 | 1.2e+03 | 0 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | TEST |  | 1.43 | 0.408 | 137 | 1.54e+03 | 1.2e+03 | 2 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | TEST |  | 6.21 | 0.0424 | 137 | 357 | 1.2e+03 | 2 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | TEST |  | 38.5 | 8e-05 | 137 | 89 | 1.2e+03 | 3 |
| baseline | B0 |  | hct116 | TEST |  | 1.57 | 0.364 | 137 | 1.4e+03 | 1.2e+03 | 2 |
| baseline | B1 |  | hct116 | TEST |  | 3.69 | 0.00661 | 137 | 1.8e+03 | 1.2e+03 | 6 |
| baseline | B2 |  | hct116 | TEST |  | 1.6 | 0.465 | 137 | 689 | 1.2e+03 | 1 |
| baseline | B3 |  | hct116 | TEST |  | 1.1e+03 | 1 | 137 | 0 | 1.2e+03 | 0 |
| baseline | B4 |  | hct116 | TEST |  | 11.2 | 4.15e-10 | 137 | 1.3e+03 | 1.2e+03 | 13 |
| baseline | B4e |  | hct116 | TEST |  | 6.6 | 4.15e-05 | 137 | 1.35e+03 | 1.2e+03 | 8 |
| baseline | B5 |  | hct116 | TEST |  | 1.31 | 1 | 137 | 419 | 1.2e+03 | 0 |
| baseline | B6 |  | hct116 | TEST |  | 1.72 | 0.324 | 137 | 1.28e+03 | 1.2e+03 | 2 |
| control | NC1 | NC1: residual-permutation null | hct116 | TEST |  | 3.35 | 0.122 | 137 | 660 | 1.2e+03 | 2 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | TEST |  | 16.9 | 4.54e-08 | 137 | 533 | 1.2e+03 | 8 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | TEST |  | 5.37 | 1 | 137 | 102 | 1.2e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | TEST |  | 6.67 | 1 | 137 | 82 | 1.2e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | TEST |  | 3.12 | 1 | 137 | 176 | 1.2e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | TEST |  | 7.1 | 1 | 137 | 77 | 1.2e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | TEST |  | 6.67 | 1 | 137 | 82 | 1.2e+03 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | paralog_exclusion | 6.66 | 1 | 137 | 82 | 1.2e+03 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_i | 5.24 | 1 | 30 | 28 | 971 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_iii | 145 | 1 | 99 | 14 | 227 | 0 |
| control | M1 |  | hct116 | TEST | NC2_M1 | 5.99 | 0.155 | 85 | 649 | 213 | 1 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | NC2 | 3.85e+03 | 1 | 85 | 0 | 213 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | no_split_sibling | 9.4 | 1 | 56 | 30 | 952 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | TEST |  | 6.75 | 1 | 137 | 81 | 1.2e+03 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | TEST |  | 15.5 | 2.36e-05 | 137 | 362 | 1.2e+03 | 5 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST |  | 4.04 | 1 | 67 | 82 | 1.05e+03 | 0 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | TEST |  | 1.09 | 0.6 | 67 | 611 | 1.05e+03 | 1 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | TEST |  | 30.5 | 6.62e-09 | 67 | 161 | 1.05e+03 | 7 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | TEST |  | 1.11 | 1 | 67 | 300 | 1.05e+03 | 0 |
| baseline | B0 |  | hek293t | TEST |  | 0.473 | 1 | 67 | 703 | 1.05e+03 | 0 |
| baseline | B1 |  | hek293t | TEST |  | 0.408 | 1 | 67 | 814 | 1.05e+03 | 0 |
| baseline | B2 |  | hek293t | TEST |  | 0.773 | 1 | 67 | 430 | 1.05e+03 | 0 |
| baseline | B3 |  | hek293t | TEST |  | 666 | 1 | 67 | 0 | 1.05e+03 | 0 |
| baseline | B4 |  | hek293t | TEST |  | 1.94 | 0.276 | 67 | 688 | 1.05e+03 | 2 |
| baseline | B4e |  | hek293t | TEST |  | 13.5 | 9.05e-12 | 67 | 714 | 1.05e+03 | 14 |
| baseline | B5 |  | hek293t | TEST |  | 283 | 4.57e-240 | 67 | 455 | 1.05e+03 | 124 |
| baseline | B6 |  | hek293t | TEST |  | 0.541 | 1 | 67 | 615 | 1.05e+03 | 0 |
| control | NC1 | NC1: residual-permutation null | hek293t | TEST |  | 0.687 | 1 | 67 | 484 | 1.05e+03 | 0 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | TEST |  | 28.3 | 4.71e-18 | 67 | 398 | 1.05e+03 | 16 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | TEST |  | 0.363 | 1 | 67 | 917 | 1.05e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | TEST |  | 4.04 | 1 | 67 | 82 | 1.05e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | TEST |  | 1.44 | 1 | 67 | 231 | 1.05e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | TEST |  | 0.297 | 1 | 67 | 1.12e+03 | 1.05e+03 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | TEST |  | 4.04 | 1 | 67 | 82 | 1.05e+03 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | paralog_exclusion | 4.03 | 1 | 67 | 82 | 1.05e+03 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_i | 23.9 | 1 | 4 | 1 | 574 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_iii | 1.35e+03 | 1 | 61 | 0 | 469 | 0 |
| control | M1 |  | hek293t | TEST | NC2_M1 | 5.15 | 1 | 45 | 216 | 210 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | NC2 | 2.23e+03 | 1 | 45 | 0 | 210 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | no_split_sibling | 4.65 | 1 | 21 | 45 | 521 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | TEST |  | 3.9 | 1 | 67 | 85 | 1.05e+03 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | TEST |  | 26.2 | 1.79e-08 | 67 | 186 | 1.05e+03 | 7 |
| control | PC2 |  | hct116 | TEST |  | 177 | 2.73e-21 | 49 | 214 | 156 | 11 |

## class:kinase

| role | run | label | line | split | subset | OR@0.1% | fisher_p(indep.) | n_regulators | selected | positives | a |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST |  | 42.9 | 1 | 51 | 69 | 83 | 0 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | TEST |  | 16.1 | 1 | 51 | 185 | 83 | 0 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | TEST |  | 36.8 | 0.0273 | 51 | 166 | 83 | 1 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | TEST |  | 61.9 | 0.0164 | 51 | 99 | 83 | 1 |
| baseline | B0 |  | hct116 | TEST |  | 5.71 | 1 | 51 | 521 | 83 | 0 |
| baseline | B1 |  | hct116 | TEST |  | 3.42 | 1 | 51 | 868 | 83 | 0 |
| baseline | B2 |  | hct116 | TEST |  | 4.05 | 1 | 51 | 734 | 83 | 0 |
| baseline | B3 |  | hct116 | TEST |  | 0.742 | 0.741 | 51 | 8.04e+03 | 83 | 1 |
| baseline | B4 |  | hct116 | TEST |  | 6.01 | 1 | 51 | 495 | 83 | 0 |
| baseline | B4e |  | hct116 | TEST |  | 5.89 | 1 | 51 | 505 | 83 | 0 |
| baseline | B5 |  | hct116 | TEST |  | 3.35 | 1 | 51 | 888 | 83 | 0 |
| baseline | B6 |  | hct116 | TEST |  | 5.36 | 1 | 51 | 555 | 83 | 0 |
| control | NC1 | NC1: residual-permutation null | hct116 | TEST |  | 4.71 | 1 | 51 | 631 | 83 | 0 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | TEST |  | 4.11 | 1 | 51 | 724 | 83 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | TEST |  | 284 | 1 | 51 | 10 | 83 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | TEST |  | 42.9 | 1 | 51 | 69 | 83 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | TEST |  | 44.8 | 1 | 51 | 66 | 83 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | TEST |  | 73.5 | 1 | 51 | 40 | 83 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | TEST |  | 39.4 | 1 | 51 | 75 | 83 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | paralog_exclusion | 43.3 | 1 | 51 | 69 | 82 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_i | 114 | 1 | 23 | 16 | 59 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_iii | 116 | 1 | 26 | 44 | 24 | 0 |
| control | M1 |  | hct116 | TEST | NC2_M1 | 20.9 | 1 | 20 | 190 | 24 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | NC2 | 7.97e+03 | 1 | 20 | 0 | 24 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | no_split_sibling | 2.11e+03 | 1 | 19 | 1 | 29 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | TEST |  | 85.4 | 0.0119 | 51 | 72 | 83 | 1 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | TEST |  | 18.8 | 1 | 51 | 158 | 83 | 0 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | TEST |  | 2.07 | 1 | 14 | 143 | 246 | 0 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| baseline | B0 |  | hek293t | TEST |  | 3.79 | 0.234 | 14 | 159 | 246 | 1 |
| baseline | B1 |  | hek293t | TEST |  | 2.27 | 0.223 | 14 | 531 | 246 | 2 |
| baseline | B2 |  | hek293t | TEST |  | 1.67 | 1 | 14 | 178 | 246 | 0 |
| baseline | B3 |  | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| baseline | B4 |  | hek293t | TEST |  | 2.03 | 1 | 14 | 146 | 246 | 0 |
| baseline | B4e |  | hek293t | TEST |  | 12.2 | 0.00217 | 14 | 151 | 246 | 3 |
| baseline | B5 |  | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| baseline | B6 |  | hek293t | TEST |  | 1.82 | 1 | 14 | 163 | 246 | 0 |
| control | NC1 | NC1: residual-permutation null | hek293t | TEST |  | 1.07 | 1 | 14 | 277 | 246 | 0 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | TEST |  | 3.33 | 1 | 14 | 89 | 246 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | paralog_exclusion | 595 | 1 | 14 | 0 | 246 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_i | 192 | 1 | 2 | 0 | 108 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_iii | 909 | 1 | 12 | 0 | 138 | 0 |
| control | M1 |  | hek293t | TEST | NC2_M1 | 3.51 | 1 | 11 | 119 | 137 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | NC2 | 840 | 1 | 11 | 0 | 137 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | no_split_sibling | 978 | 1 | 7 | 0 | 75 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | TEST |  | 596 | 1 | 14 | 0 | 246 | 0 |
| control | PC2 |  | hct116 | TEST |  | 469 | 0.00261 | 12 | 51 | 6 | 1 |

## class:other

| role | run | label | line | split | subset | OR@0.1% | fisher_p(indep.) | n_regulators | selected | positives | a |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST |  | 3.21 | 8.09e-23 | 1.6e+03 | 16429 | 29978 | 101 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | TEST |  | 1.05 | 0.423 | 1.6e+03 | 15353 | 29978 | 31 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | TEST |  | 17.2 | 0 | 1.6e+03 | 16489 | 29978 | 523 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | TEST |  | 4.01 | 1.76e-36 | 1.6e+03 | 16323 | 29978 | 125 |
| baseline | B0 |  | hct116 | TEST |  | 1.14 | 0.243 | 1.6e+03 | 15450 | 29978 | 34 |
| baseline | B1 |  | hct116 | TEST |  | 1.16 | 0.219 | 1.6e+03 | 14747 | 29978 | 33 |
| baseline | B2 |  | hct116 | TEST |  | 5.53 | 1.79e-64 | 1.6e+03 | 15387 | 29978 | 162 |
| baseline | B3 |  | hct116 | TEST |  | 0.0268 | 1 | 1.6e+03 | 9.65e+03 | 29978 | 0 |
| baseline | B4 |  | hct116 | TEST |  | 5.87 | 8.31e-73 | 1.6e+03 | 15582 | 29978 | 174 |
| baseline | B4e |  | hct116 | TEST |  | 7.74 | 2.98e-117 | 1.6e+03 | 15517 | 29978 | 227 |
| baseline | B5 |  | hct116 | TEST |  | 27 | 0 | 1.6e+03 | 15182 | 29978 | 736 |
| baseline | B6 |  | hct116 | TEST |  | 0.566 | 0.996 | 1.6e+03 | 15563 | 29978 | 17 |
| control | NC1 | NC1: residual-permutation null | hct116 | TEST |  | 1.31 | 0.0527 | 1.6e+03 | 16189 | 29978 | 41 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | TEST |  | 11.1 | 6.57e-213 | 1.6e+03 | 15746 | 29978 | 327 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | TEST |  | 2.63 | 3.7e-14 | 1.6e+03 | 16271 | 29978 | 82 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | TEST |  | 3.11 | 2.65e-21 | 1.6e+03 | 16431 | 29978 | 98 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | TEST |  | 2.03 | 3.35e-07 | 1.6e+03 | 16155 | 29978 | 63 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | TEST |  | 2.82 | 6.16e-17 | 1.6e+03 | 16479 | 29978 | 89 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | TEST |  | 3.08 | 7.87e-21 | 1.6e+03 | 16418 | 29978 | 97 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | paralog_exclusion | 3.18 | 2.48e-22 | 1.6e+03 | 16421 | 29961 | 100 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_i | 1.13 | 0.279 | 645 | 6.49e+03 | 25742 | 30 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_iii | 0.212 | 0.991 | 898 | 9.98e+03 | 4.14e+03 | 1 |
| control | M1 |  | hct116 | TEST | NC2_M1 | 1.26 | 0.368 | 781 | 7.82e+03 | 3.88e+03 | 5 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | NC2 | 0.113 | 1 | 781 | 8.7e+03 | 3.88e+03 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | no_split_sibling | 0.894 | 0.674 | 637 | 6.36e+03 | 9.89e+03 | 9 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | TEST |  | 2.54 | 3.85e-13 | 1.6e+03 | 16407 | 29978 | 80 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | TEST |  | 16.4 | 0 | 1.6e+03 | 16472 | 29978 | 498 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST |  | 0.509 | 0.977 | 710 | 8.28e+03 | 10604 | 6 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | TEST |  | 0.469 | 0.981 | 710 | 7.49e+03 | 10604 | 5 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | TEST |  | 4.78 | 5.83e-20 | 710 | 8.01e+03 | 10604 | 54 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | TEST |  | 0.0869 | 1 | 710 | 8.07e+03 | 10604 | 1 |
| baseline | B0 |  | hek293t | TEST |  | 0.66 | 0.904 | 710 | 7.45e+03 | 10604 | 7 |
| baseline | B1 |  | hek293t | TEST |  | 1 | 0.54 | 710 | 7.02e+03 | 10604 | 10 |
| baseline | B2 |  | hek293t | TEST |  | 1.55 | 0.0556 | 710 | 7.73e+03 | 10604 | 17 |
| baseline | B3 |  | hek293t | TEST |  | 0.755 | 0.839 | 710 | 8.37e+03 | 10604 | 9 |
| baseline | B4 |  | hek293t | TEST |  | 0.941 | 0.618 | 710 | 7.47e+03 | 10604 | 10 |
| baseline | B4e |  | hek293t | TEST |  | 14.7 | 1.61e-116 | 710 | 7.44e+03 | 10604 | 151 |
| baseline | B5 |  | hek293t | TEST |  | 29.5 | 9.95e-292 | 710 | 7.08e+03 | 10604 | 279 |
| baseline | B6 |  | hek293t | TEST |  | 0.47 | 0.981 | 710 | 7.47e+03 | 10604 | 5 |
| control | NC1 | NC1: residual-permutation null | hek293t | TEST |  | 1.95 | 0.00381 | 710 | 7.6e+03 | 10604 | 21 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | TEST |  | 9.08 | 1.46e-56 | 710 | 7.58e+03 | 10604 | 96 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | TEST |  | 1.7 | 0.0235 | 710 | 7.45e+03 | 10604 | 18 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | TEST |  | 0.509 | 0.977 | 710 | 8.28e+03 | 10604 | 6 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | TEST |  | 0.611 | 0.938 | 710 | 8.04e+03 | 10604 | 7 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | TEST |  | 0.582 | 0.944 | 710 | 7.24e+03 | 10604 | 6 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | TEST |  | 0.509 | 0.977 | 710 | 8.28e+03 | 10604 | 6 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | paralog_exclusion | 0.509 | 0.977 | 710 | 8.28e+03 | 10603 | 6 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_i | 0.559 | 0.872 | 62 | 754 | 3.08e+03 | 2 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_iii | 0.0595 | 1 | 641 | 7.52e+03 | 7.5e+03 | 0 |
| control | M1 |  | hek293t | TEST | NC2_M1 | 0.329 | 0.984 | 544 | 5.96e+03 | 5.82e+03 | 2 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | NC2 | 0.0773 | 1 | 544 | 6.33e+03 | 5.82e+03 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | no_split_sibling | 0.121 | 1 | 288 | 3.32e+03 | 3.79e+03 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | TEST |  | 0.764 | 0.83 | 710 | 8.28e+03 | 10604 | 9 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | TEST |  | 4.43 | 2.55e-17 | 710 | 7.99e+03 | 10604 | 50 |
| control | PC2 |  | hct116 | TEST |  | 177 | 0 | 531 | 5.36e+03 | 1.87e+03 | 277 |

## class:ribosomal

| role | run | label | line | split | subset | OR@0.1% | fisher_p(indep.) | n_regulators | selected | positives | a |
|---|---|---|---|---|---|---|---|---|---|---|---|
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST |  | 5.43 | 1 | 16 | 631 | 22 | 0 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hct116 | TEST |  | 8.69 | 1 | 16 | 395 | 22 | 0 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hct116 | TEST |  | 8.17 | 1 | 16 | 420 | 22 | 0 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hct116 | TEST |  | 3.75 | 1 | 16 | 913 | 22 | 0 |
| baseline | B0 |  | hct116 | TEST |  | 20.3 | 1 | 16 | 169 | 22 | 0 |
| baseline | B1 |  | hct116 | TEST |  | 45.6 | 1 | 16 | 75 | 22 | 0 |
| baseline | B2 |  | hct116 | TEST |  | 147 | 1 | 16 | 23 | 22 | 0 |
| baseline | B3 |  | hct116 | TEST |  | 6.89e+03 | 1 | 16 | 0 | 22 | 0 |
| baseline | B4 |  | hct116 | TEST |  | 22.4 | 1 | 16 | 153 | 22 | 0 |
| baseline | B4e |  | hct116 | TEST |  | 22.1 | 1 | 16 | 155 | 22 | 0 |
| baseline | B5 |  | hct116 | TEST |  | 3.36 | 1 | 16 | 1.02e+03 | 22 | 0 |
| baseline | B6 |  | hct116 | TEST |  | 22.9 | 1 | 16 | 150 | 22 | 0 |
| control | NC1 | NC1: residual-permutation null | hct116 | TEST |  | 23.2 | 1 | 16 | 148 | 22 | 0 |
| control | C2_g2 | C2: context fusion applied to NC1 | hct116 | TEST |  | 6.62 | 1 | 16 | 518 | 22 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hct116 | TEST |  | 5.05 | 1 | 16 | 679 | 22 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hct116 | TEST |  | 5.5 | 1 | 16 | 624 | 22 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hct116 | TEST |  | 5.03 | 1 | 16 | 682 | 22 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hct116 | TEST |  | 5.36 | 1 | 16 | 640 | 22 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hct116 | TEST |  | 5.43 | 1 | 16 | 631 | 22 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | paralog_exclusion | 5.43 | 1 | 16 | 631 | 22 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_i | 11.3 | 1 | 15 | 285 | 22 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | class_iii | 19747 | 1 | 1 | 0 | 0 | 0 |
| control | M1 |  | hct116 | TEST | NC2_M1 | 858 | 1 | 1 | 11 | 0 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | NC2 | 19747 | 1 | 1 | 0 | 0 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hct116 | TEST | no_split_sibling | 6.23 | 1 | 6 | 347 | 13 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hct116 | TEST |  | 5.37 | 1 | 16 | 639 | 22 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hct116 | TEST |  | 7.68 | 1 | 16 | 447 | 22 | 0 |
| selected | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST |  | 362 | 1 | 3 | 12 | 3 | 0 |
| ablation | a3_r095_s2_w0 | ablation: selected minus G1 | hek293t | TEST |  | 36.8 | 1 | 3 | 122 | 3 | 0 |
| ablation | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2 | context layer (selected + I1 γ=2) | hek293t | TEST |  | 22.1 | 1 | 3 | 203 | 3 | 0 |
| ablation | focal__g1_t100_w0__a3_r095_s2_w0 | ablation: focal arm handling | hek293t | TEST |  | 1.29e+03 | 1 | 3 | 3 | 3 | 0 |
| baseline | B0 |  | hek293t | TEST |  | 153 | 1 | 3 | 29 | 3 | 0 |
| baseline | B1 |  | hek293t | TEST |  | 335 | 1 | 3 | 13 | 3 | 0 |
| baseline | B2 |  | hek293t | TEST |  | 159 | 1 | 3 | 28 | 3 | 0 |
| baseline | B3 |  | hek293t | TEST |  | 9.06e+03 | 1 | 3 | 0 | 3 | 0 |
| baseline | B4 |  | hek293t | TEST |  | 139 | 1 | 3 | 32 | 3 | 0 |
| baseline | B4e |  | hek293t | TEST |  | 139 | 1 | 3 | 32 | 3 | 0 |
| baseline | B5 |  | hek293t | TEST |  | 5.25 | 1 | 3 | 839 | 3 | 0 |
| baseline | B6 |  | hek293t | TEST |  | 55.4 | 1 | 3 | 81 | 3 | 0 |
| control | NC1 | NC1: residual-permutation null | hek293t | TEST |  | 697 | 1 | 3 | 6 | 3 | 0 |
| control | C2_g2 | C2: context fusion applied to NC1 | hek293t | TEST |  | 24.1 | 1 | 3 | 186 | 3 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__S1 |  | hek293t | TEST |  | 697 | 1 | 3 | 6 | 3 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__noHCT116 |  | hek293t | TEST |  | 362 | 1 | 3 | 12 | 3 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__hidden_factors |  | hek293t | TEST |  | 312 | 1 | 3 | 14 | 3 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__cn_outlier_exclusion |  | hek293t | TEST |  | 431 | 1 | 3 | 10 | 3 | 0 |
| sensitivity | g1_t100_w0__a3_r095_s2_w0__tp53_chr17_drop |  | hek293t | TEST |  | 362 | 1 | 3 | 12 | 3 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | paralog_exclusion | 362 | 1 | 3 | 12 | 3 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_i | 6.03e+03 | 1 | 2 | 0 | 3 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | class_iii | 21233 | 1 | 1 | 0 | 0 | 0 |
| control | M1 |  | hek293t | TEST | NC2_M1 | 1.12e+03 | 1 | 1 | 9 | 0 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | NC2 | 21233 | 1 | 1 | 0 | 0 | 0 |
| subset | g1_t100_w0__a3_r095_s2_w0 | selected (CN-only: M1 + G1 + A3) | hek293t | TEST | no_split_sibling | 824 | 1 | 3 | 5 | 3 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0_C1 |  | hek293t | TEST |  | 394 | 1 | 3 | 11 | 3 | 0 |
| control | g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1 |  | hek293t | TEST |  | 23 | 1 | 3 | 195 | 3 | 0 |
| control | PC2 |  | hct116 | TEST |  | 8.42e+03 | 1 | 3 | 0 | 3 | 0 |
