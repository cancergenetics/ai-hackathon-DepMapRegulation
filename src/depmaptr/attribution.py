"""Co-segment attribution A3 (DESIGN §5.4).

1. Greedy sequential segmentation along each chromosome: a gene joins the current segment iff
   corr(x_g, mean profile of the current segment) >= rho; otherwise it starts a new segment.
   x is the active-transform CN residualised on the LOCO design D_c (implementation choice, docs/MODEL.md).
2. Within a segment (ALL CN genes in it), single-effect (SuSiE L = 1) posterior per target B:
       alpha_{g|B} ∝ pi_g * ABF(z_gB) * ABF(z_co,gB)^w_co
   ABF(z) = sqrt(1 - lam) exp(lam z^2 / 2), lam = sigma0^2 / (sigma0^2 + 1/df), sigma0 = 0.05;
   pi_g ∝ exp(t_cis,g / s) * w_G2(g); if sum(pi) = 0 within a segment, pi is uniform.
   Singleton segments get alpha = 1. Computed in log space.
3. z_co,gB = atanh(r) sqrt(df - 1), r = partial corr(y_g, y_B | x_g, D, W_B) (computed in pipeline).
"""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

SIGMA0 = 0.05


def greedy_segments(xres: np.ndarray, rho: float) -> np.ndarray:
    """xres: n x m CN (columns in genomic order). Returns segment labels 0..S-1 (length m)."""
    n, m = xres.shape
    labels = np.empty(m, np.int64)
    if m == 0:
        return labels
    X = xres - xres.mean(0)
    seg = 0
    labels[0] = 0
    ssum = X[:, 0].copy()
    for j in range(1, m):
        mu = ssum  # mean profile up to a positive scale; correlation is scale-invariant
        a = X[:, j]
        den = np.sqrt((a @ a) * (mu @ mu))
        r = (a @ mu) / den if den > 0 else 0.0
        if r >= rho:
            labels[j] = seg
            ssum += a
        else:
            seg += 1
            labels[j] = seg
            ssum = a.copy()
    return labels


def log_abf(z: np.ndarray, df, sigma0: float = SIGMA0) -> np.ndarray:
    lam = sigma0 ** 2 / (sigma0 ** 2 + 1.0 / np.asarray(df, np.float64))
    return 0.5 * np.log1p(-lam) + 0.5 * lam * np.asarray(z, np.float64) ** 2


def segment_alpha(labels: np.ndarray, t: np.ndarray, df: np.ndarray, t_cis: np.ndarray, s: float,
                  w_g2: np.ndarray | None = None, zco: np.ndarray | None = None, w_co: float = 0.0,
                  df_co: np.ndarray | None = None, sigma0: float = SIGMA0) -> np.ndarray:
    """alpha (m x T) for all genes of one chromosome (rows in the order of `labels`).

    t, df: m x T reduced-form t and df of each segment gene vs each target (NaN t -> z = 0).
    t_cis: m; w_g2: m (None -> 1); zco, df_co: m x T (only used when w_co > 0).
    """
    m, T = t.shape
    z = np.nan_to_num(t, nan=0.0)
    ll = log_abf(z, df, sigma0)
    if w_co > 0:
        assert zco is not None and df_co is not None
        ll = ll + w_co * log_abf(np.nan_to_num(zco, nan=0.0), df_co, sigma0)
    with np.errstate(divide="ignore"):
        lp = np.nan_to_num(np.asarray(t_cis, np.float64), nan=0.0) / s
        if w_g2 is not None:
            lp = lp + np.log(np.asarray(w_g2, np.float64))
    alpha = np.ones((m, T), np.float64)
    for seg in np.unique(labels):
        idx = np.flatnonzero(labels == seg)
        if len(idx) == 1:
            continue
        lpi = lp[idx]
        if not np.isfinite(lpi).any():          # sum(pi) = 0 -> uniform prior
            lpi = np.zeros(len(idx))
        a = lpi[:, None] + ll[idx]
        alpha[idx] = np.exp(a - logsumexp(a, axis=0, keepdims=True))
    return alpha
