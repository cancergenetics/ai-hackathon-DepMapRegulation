"""Generic chunked float64 association engine (DESIGN §5.1; algorithms review #1-#3, #11).

Model for a pair (A, B) after Frisch-Waugh-Lovell on a shared design D (residualised by pivoted QR):

    y_B ~ u_A + V_A + W_B          (all columns residualised on D)

- u_A: the tested regulator-specific regressor (x_A; the loss part of the hinge; y_g for co-expression);
- V_A: regulator-specific nuisance regressors (gain part of the hinge; x_g for co-expression);
- W_B: target-specific nuisance regressors (cis basis [x_B, 1(relCN_B < 0.25)]); an all-zero column means
  "absent for this target" and is skipped, so k_B and df are per target.

For every pair the (k+2)x(k+2) Gram matrix of [V, W, u, y] is assembled from GEMMs (regulator x target
cross-products) and broadcast vectors (regulator-only / target-only products), then swept elementwise
(batched Gaussian elimination). The partial correlation of u and y is S_uy / sqrt(S_uu S_yy) and

    t = r sqrt(df / (1 - r^2)),   df = n - rank(D) - (#non-degenerate nuisance pivots) - 1.

This equals the OLS t of u in y ~ [D, V, W, u] (tested against statsmodels to 1e-8).
Guards (DESIGN §5.1): NaN where the residual norm of u or y (after D, V, W) falls below 1e-3 of its
centred raw norm, or df < MIN_DF.
"""
from __future__ import annotations

import numpy as np
from scipy import linalg

NORM_GUARD = 1e-3
MIN_DF = 3
PIVOT_TOL = 1e-10
COND_MAX = 1e6
RANK_TOL = 1e-9


class Residualiser:
    """Pivoted-QR projection off span(D) in float64; logs rank and condition number."""

    def __init__(self, D: np.ndarray, cond_max: float = COND_MAX):
        D = np.asarray(D, np.float64)
        sc = np.linalg.norm(D, axis=0)
        keep = sc > 0
        Ds = D[:, keep] / sc[keep]
        Q, R, P = linalg.qr(Ds, mode="economic", pivoting=True)
        d = np.abs(np.diag(R))
        rank = int((d > RANK_TOL * d[0]).sum())
        self.rank = rank
        self.n_cols = D.shape[1]
        self.cond = float(np.linalg.cond(R[:rank, :rank]))
        if not np.isfinite(self.cond) or self.cond >= cond_max:
            raise AssertionError(f"cond(D) = {self.cond:.3g} >= {cond_max:g} (rank {rank}/{D.shape[1]})")
        self.Q = np.ascontiguousarray(Q[:, :rank])
        self.n = D.shape[0]

    def __call__(self, M: np.ndarray) -> np.ndarray:
        M = np.asarray(M, np.float64)
        return M - self.Q @ (self.Q.T @ M)

    @property
    def df_base(self) -> int:
        return self.n - self.rank


def centred_norm(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, np.float64)
    return np.linalg.norm(M - M.mean(0), axis=0)


def t_from_r(r: np.ndarray, df) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return r * np.sqrt(df / np.maximum(1.0 - r * r, 1e-300))


def pair_partial(U: np.ndarray, Y: np.ndarray, V=(), W=(), *, df_base: int,
                 u_raw: np.ndarray | None = None, y_raw: np.ndarray | None = None,
                 norm_guard: float = NORM_GUARD, min_df: int = MIN_DF, pivot_tol: float = PIVOT_TOL):
    """Partial correlation r, t and df of u_A with y_B given V_A and W_B, for all (A, B).

    U: n x R, Y: n x T, V: list of n x R, W: list of n x T, all already residualised on D.
    u_raw (R,) / y_raw (T,): centred raw norms for the residual-norm guard (None: no guard).
    Returns r, t, df as float64 arrays of shape (R, T).
    """
    V, W = list(V), list(W)
    cols = [("R", v) for v in V] + [("T", w) for w in W] + [("R", U), ("T", Y)]
    k = len(cols)

    def gram(i, j):
        (ki, a), (kj, b) = cols[i], cols[j]
        if ki == "R" and kj == "R":
            return np.einsum("ij,ij->j", a, b)[:, None]
        if ki == "T" and kj == "T":
            return np.einsum("ij,ij->j", a, b)[None, :]
        if ki == "R":
            return a.T @ b          # (R, T)
        return b.T @ a              # (R, T): all mixed blocks are oriented regulator x target

    G = {(i, j): gram(i, j) for i in range(k) for j in range(i, k)}
    npiv = 0
    for p in range(k - 2):
        orig = G[(p, p)]
        d = orig
        # orig diag before this sweep vs. after previous sweeps: degenerate if collapsed
        base = np.einsum("ij,ij->j", cols[p][1], cols[p][1])
        base = base[:, None] if cols[p][0] == "R" else base[None, :]
        ok = (base > 0) & (d > pivot_tol * base)
        with np.errstate(divide="ignore", invalid="ignore"):
            inv = np.where(ok, 1.0 / np.where(ok, d, 1.0), 0.0)
        npiv = npiv + ok
        for i in range(p + 1, k):
            gip = G[(p, i)]
            f = gip * inv
            for j in range(i, k):
                G[(i, j)] = G[(i, j)] - f * G[(p, j)]
    Suu, Suy, Syy = G[(k - 2, k - 2)], G[(k - 2, k - 1)], G[(k - 1, k - 1)]
    Suu = np.maximum(Suu, 0.0)
    Syy = np.maximum(Syy, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = Suy / np.sqrt(Suu * Syy)
    df = df_base - npiv - 1
    r, df = np.broadcast_arrays(r, df)
    r = np.array(r, np.float64)
    df = np.array(df, np.float64)
    bad = ~np.isfinite(r) | (df < min_df)
    if u_raw is not None:
        bad |= np.broadcast_to(np.sqrt(Suu) < norm_guard * np.asarray(u_raw)[:, None], r.shape)
    if y_raw is not None:
        bad |= np.broadcast_to(np.sqrt(Syy) < norm_guard * np.asarray(y_raw)[None, :], r.shape)
    r[bad] = np.nan
    r = np.clip(r, -1.0, 1.0)
    t = t_from_r(r, df)
    return r, t, df


def diag_partial(U: np.ndarray, Y: np.ndarray, *, df_base: int, u_raw=None, y_raw=None,
                 norm_guard: float = NORM_GUARD, min_df: int = MIN_DF):
    """Column-wise partial correlation of U[:, j] with Y[:, j] (both residualised on D); df = n - rank - 1."""
    Suu = np.einsum("ij,ij->j", U, U)
    Syy = np.einsum("ij,ij->j", Y, Y)
    Suy = np.einsum("ij,ij->j", U, Y)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = Suy / np.sqrt(Suu * Syy)
    df = float(df_base - 1)
    bad = ~np.isfinite(r) | (df < min_df)
    if u_raw is not None:
        bad |= np.sqrt(Suu) < norm_guard * np.asarray(u_raw)
    if y_raw is not None:
        bad |= np.sqrt(Syy) < norm_guard * np.asarray(y_raw)
    r = np.where(bad, np.nan, r)
    return r, t_from_r(r, df), df


def chunked_pair_partial(U, Y, V=(), W=(), *, df_base, u_raw=None, y_raw=None, chunk=512, **kw):
    """pair_partial over regulator chunks (bounds peak memory: ~ (k+2)^2/2 R_chunk x T float64 arrays)."""
    R = U.shape[1]
    T = Y.shape[1]
    r_out = np.empty((R, T), np.float64)
    t_out = np.empty((R, T), np.float64)
    df_out = np.empty((R, T), np.float64)
    for s in range(0, R, chunk):
        e = min(R, s + chunk)
        r, t, df = pair_partial(U[:, s:e], Y, [v[:, s:e] for v in V], W, df_base=df_base,
                                u_raw=None if u_raw is None else u_raw[s:e], y_raw=y_raw, **kw)
        r_out[s:e], t_out[s:e], df_out[s:e] = r, t, df
    return r_out, t_out, df_out


def lineage_permutation(strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Row permutation within strata (NC1, DESIGN §5.2): perm[i] is the source row for row i."""
    strata = np.asarray(strata)
    perm = np.arange(len(strata))
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s)
        perm[idx] = rng.permutation(idx)
    return perm


def resid_ss(Y: np.ndarray, W=(), pivot_tol: float = PIVOT_TOL) -> np.ndarray:
    """Per column j: residual sum of squares of Y[:, j] after regressing on the W_k[:, j] (all residualised on D).
    All-zero (inactive) W columns are skipped, as in pair_partial."""
    cols = list(W) + [Y]
    k = len(cols)
    G = {(i, j): np.einsum("ij,ij->j", cols[i], cols[j]) for i in range(k) for j in range(i, k)}
    for p in range(k - 1):
        base = np.einsum("ij,ij->j", cols[p], cols[p])
        d = G[(p, p)]
        ok = (base > 0) & (d > pivot_tol * base)
        inv = np.where(ok, 1.0 / np.where(ok, d, 1.0), 0.0)
        for i in range(p + 1, k):
            f = G[(p, i)] * inv
            for j in range(i, k):
                G[(i, j)] = G[(i, j)] - f * G[(p, j)]
    return np.maximum(G[(k - 1, k - 1)], 0.0)
