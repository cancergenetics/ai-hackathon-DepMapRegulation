"""Gates (DESIGN §1, §5.3).

G1 (instrument): the cis t of y_A ~ x_A | D_g, where D_g is the design of A's group under the active arm
handling (computed in pipeline.py with assoc.diag_partial). A passes iff t_cis > 0 and t_cis^2 >= tau.
Regulators without an expression column have t_cis = 0 (fail G1, s_A = +1).

s_A (knockdown sign rule, DESIGN §1): s_A = -1 iff t_cis <= -3, else +1.
Predicted knockdown sign of B: -sign(beta_AB) * s_A (harness; `predicted_kd_sign` here for tests).

G2 (context expression): A is expressed in the context line (DepMap ACH-000971 for HCT116,
ACH-001085 for HEK293T; IDs asserted against Model.csv in io.context_model_check) iff
log2(TPM+1) >= threshold. Regulators without an expression column count as not expressed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

S_CIS_THRESHOLD = -3.0


def s_cis(t_cis: np.ndarray) -> np.ndarray:
    t = np.nan_to_num(np.asarray(t_cis, np.float64), nan=0.0)
    return np.where(t <= S_CIS_THRESHOLD, -1, 1).astype(np.int8)


def predicted_kd_sign(score: np.ndarray, s: np.ndarray) -> np.ndarray:
    """-sign(score) * s_A, broadcasting s over targets."""
    return -np.sign(score) * np.asarray(s)[:, None]


def g1_pass(t_cis: np.ndarray, tau: float) -> np.ndarray:
    t = np.nan_to_num(np.asarray(t_cis, np.float64), nan=0.0)
    return (t > 0) & (t * t >= tau)


def g1_weight(t_cis: np.ndarray, tau: float, w: float) -> np.ndarray:
    return np.where(g1_pass(t_cis, tau), 1.0, w)


def context_expressed(regulators, context_expr: pd.Series, threshold: float) -> np.ndarray:
    v = context_expr.reindex(list(regulators)).to_numpy(np.float64)
    return np.nan_to_num(v, nan=-np.inf) >= threshold


def g2_weight(regulators, context_expr: pd.Series, threshold: float, w: float) -> np.ndarray:
    return np.where(context_expressed(regulators, context_expr, threshold), 1.0, w)
