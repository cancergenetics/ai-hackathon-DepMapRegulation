"""Baselines B0-B6 (DESIGN §7). Each returns the same dict shape as pipeline.run (regulators x targets).

Sign conventions (s_cis): B2 is an mRNA-level statistic, so s_cis = +1 for every regulator (lowering A's
mRNA predicts -sign(r)). B1 uses its own naive cis t (no covariates). B0, B3-B6 use the M1 (LOCO) G1 cis t
for s_cis; for the target-agnostic baselines the sign carries no information.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata

from . import gates
from .assoc import t_from_r
from .io import Dataset, cn_transform
from .pipeline import Params, cis_t, prepare, run_arm_only

ROOT_SEED = 20261027


def _std(M):
    M = np.asarray(M, np.float64)
    M = M - M.mean(0)
    nrm = np.linalg.norm(M, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return M / nrm


def _pearson_t(A, B, chunk=1024):
    """Pearson r -> t (df = n - 2) between all columns of A (n x R) and B (n x T)."""
    n = A.shape[0]
    As, Bs = _std(A), _std(B)
    out = np.empty((A.shape[1], B.shape[1]), np.float32)
    for s in range(0, A.shape[1], chunk):
        r = np.clip(As[:, s:s + chunk].T @ Bs, -1, 1)
        out[s:s + chunk] = t_from_r(r, n - 2)
    return out


def _idx(keys, ref):
    pos = {k: i for i, k in enumerate(ref)}
    return np.array([pos.get(k, -1) for k in keys])


def m1_cis(ds: Dataset, p: Params | None = None) -> np.ndarray:
    p = p or Params()
    return cis_t(prepare(ds, p), np.asarray(ds.regulators))


def _pack(t, cis, s=None, diag=None):
    t = np.asarray(t, np.float32)
    cis = np.asarray(cis, np.float64)
    return {"t": t, "score": t.copy(), "cis_t": cis.astype(np.float32),
            "s_cis": gates.s_cis(cis) if s is None else s, "diagnostics": diag or {}}


def run_baseline(ds: Dataset, bid: str, p: Params | None = None, seed: int | None = None, log=print) -> dict:
    p = p or Params()
    R, T = len(ds.regulators), len(ds.targets)
    rc = _idx(ds.regulators, ds.cn_keys)
    tc = _idx(ds.targets, ds.cn_keys)
    re_ = _idx(ds.regulators, ds.expr_keys)
    te = _idx(ds.targets, ds.expr_keys)
    x = cn_transform(ds.relcn, floor=p.cn_floor, cap=p.cn_cap)
    if bid == "B0":
        ss = np.random.SeedSequence(ROOT_SEED if seed is None else seed)
        t = np.random.default_rng(ss).random((R, T), dtype=np.float32)
        return _pack(t, m1_cis(ds, p), diag={"seed": int(ROOT_SEED if seed is None else seed)})
    if bid == "B1":
        t = _pearson_t(x[:, rc], ds.expr[:, te])
        # naive cis t (no covariates) for the sign rule
        cis = np.zeros(R)
        h = re_ >= 0
        xa, ya = _std(x[:, rc[h]]), _std(ds.expr[:, re_[h]])
        cis[h] = np.nan_to_num(t_from_r(np.clip((xa * ya).sum(0), -1, 1), ds.n - 2), nan=0.0)
        return _pack(t, cis)
    if bid == "B2":
        t = np.full((R, T), np.nan, np.float32)
        h = re_ >= 0
        t[h] = _pearson_t(ds.expr[:, re_[h]], ds.expr[:, te])
        # Same-gene entries (A == B) are r = 1 -> t = inf; they are masked by the harness (same chromosome).
        return _pack(t, m1_cis(ds, p), s=np.ones(R, np.int8), diag={"n_regulators_without_expression": int((~h).sum())})
    cis = m1_cis(ds, p)
    if bid == "B3":
        sd = x[:, rc].std(0, ddof=1)
        return _pack(np.broadcast_to((cis * sd)[:, None], (R, T)).astype(np.float32), cis)
    if bid == "B4":
        sd = ds.expr[:, te].std(0, ddof=1)
        return _pack(np.broadcast_to(sd[None, :], (R, T)).astype(np.float32), cis)
    if bid in ("B4e", "B5"):
        ctx = ds.context_expr[p.context].reindex(list(ds.targets)).to_numpy(np.float64)
        assert np.isfinite(ctx).all(), "context expression missing for some targets"
        if bid == "B4e":
            return _pack(np.broadcast_to(ctx[None, :], (R, T)).astype(np.float32), cis)
        sd = x[:, rc].std(0, ddof=1)
        rank_a = rankdata(np.abs(cis * sd))          # B3 ranks by |score|
        rank_b = rankdata(ctx)
        return _pack(np.outer(rank_a, rank_b).astype(np.float32), cis)
    if bid == "B6":
        return run_arm_only(ds, p, log=log)
    raise ValueError(f"unknown baseline {bid}")
