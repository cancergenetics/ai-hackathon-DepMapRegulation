"""Calibration diagnostics (DESIGN §5.2): genomic inflation lambda_GC from t statistics."""
from __future__ import annotations

import numpy as np
from scipy import special, stats

CHI2_MEDIAN_1DF = float(stats.chi2.ppf(0.5, 1))  # 0.4549...


def z_from_t(t: np.ndarray, df) -> np.ndarray:
    """Signed z with the same two-sided tail probability as t (log-space; never underflows to |z| = inf)."""
    t = np.asarray(t, np.float64)
    df = np.broadcast_to(np.asarray(df, np.float64), t.shape)
    logp = stats.t.logsf(np.abs(t), df)            # one-sided
    with np.errstate(divide="ignore", invalid="ignore"):
        z = -special.ndtri_exp(logp)
        # where the t tail underflows, use Wallace's (1959) tail approximation z = sqrt(df log(1 + t^2/df))
        wal = np.sqrt(df * np.log1p(t * t / df)) * (1 - 1 / (8 * df))
    z = np.where(np.isfinite(z), z, wal)
    return np.sign(t) * z


def lambda_gc(t: np.ndarray, df=None) -> float:
    t = np.asarray(t, np.float64).ravel()
    ok = np.isfinite(t)
    if df is None:
        z = t[ok]
    else:
        z = z_from_t(t[ok], np.broadcast_to(np.asarray(df, np.float64), np.asarray(t).shape).ravel()[ok])
    return float(np.median(z * z) / CHI2_MEDIAN_1DF)
