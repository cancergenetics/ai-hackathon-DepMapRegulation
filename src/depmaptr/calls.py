"""MVP call set (DESIGN §12, until the Efron empirical null exists):

    call = {trans pairs with |t| >= t*},
    t* = the smallest threshold with (#NC1 trans pairs with |t| >= t* / n_perm) / (#real trans pairs with |t| >= t*) <= 0.05.

The NC1 passes must come from the same variant (same t-affecting parameters). Counts use a fine histogram
of |t| (bin width 1e-3 up to 100; values above go to the top bin), so t* is a bin edge. `n_calls` in the
score-file metadata is recounted from the call matrix itself.
"""
from __future__ import annotations

import numpy as np

BIN = 1e-3
NBINS = 100_000
EDGES = np.arange(NBINS) * BIN                          # lower edges; the last bin is open-ended
FDR = 0.05


def trans_mask(reg_chrom: np.ndarray, tgt_chrom: np.ndarray) -> np.ndarray:
    return np.asarray(reg_chrom)[:, None] != np.asarray(tgt_chrom)[None, :]


def abs_t_hist(t: np.ndarray, trans: np.ndarray) -> np.ndarray:
    """Histogram of |t| over finite trans entries, processed in row blocks."""
    h = np.zeros(NBINS, np.int64)
    for s in range(0, t.shape[0], 1024):
        blk = np.abs(t[s:s + 1024][trans[s:s + 1024]]).astype(np.float64)
        blk = blk[np.isfinite(blk)]
        idx = np.minimum((blk / BIN).astype(np.int64), NBINS - 1)
        h += np.bincount(idx, minlength=NBINS)
    return h


def t_star(real_hist: np.ndarray, null_hist: np.ndarray, n_perm: int, fdr: float = FDR) -> dict:
    real_ge = np.cumsum(real_hist[::-1])[::-1]          # #|t| >= EDGES[k]
    null_ge = np.cumsum(null_hist[::-1])[::-1] / n_perm
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(real_ge > 0, null_ge / real_ge, np.inf)
    ok = np.flatnonzero(ratio <= fdr)
    if len(ok) == 0:
        return {"t_star": float("inf"), "n_calls": 0, "fdr_at_t_star": None, "n_perm": int(n_perm)}
    k = int(ok[0])
    return {"t_star": float(EDGES[k]), "n_calls": int(real_ge[k]), "fdr_at_t_star": float(ratio[k]),
            "n_perm": int(n_perm), "n_null_ge_t_star": float(null_ge[k])}


def call_matrix(t: np.ndarray, trans: np.ndarray, tstar: float) -> np.ndarray:
    with np.errstate(invalid="ignore"):
        return trans & (np.abs(t) >= tstar)
