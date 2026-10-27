"""Known-answer tests for the association engine (DESIGN §5.1; algorithms review #3):
closed-form FWL partial-correlation t vs brute-force statsmodels OLS to 1e-8, on toy data and on >= 200
random real pairs (production code path), including the homozygous-deletion indicator, per-target df,
the hinge (regulator-specific nuisance), the co-expression form and a rank-deficient design."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm

from depmaptr.assoc import Residualiser, chunked_pair_partial, diag_partial, pair_partial
from depmaptr.calibrate import lambda_gc, z_from_t

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
TOL = 1e-8


def _ols_t(y, cols):
    f = sm.OLS(y, np.column_stack(cols)).fit()
    return f.tvalues[-1], f.df_resid


def _toy(seed=0, n=300, p=10, R=12, T=20, rank_deficient=True):
    rng = np.random.default_rng(seed)
    D = np.column_stack([np.ones(n), rng.normal(size=(n, p)), rng.integers(0, 2, (n, 4))])
    if rank_deficient:
        D = np.column_stack([D, D[:, 1] - 2 * D[:, 3]])
    X = rng.normal(size=(n, R))
    XB = rng.normal(size=(n, T)) + 0.3 * X[:, :1]
    DEL = np.zeros((n, T))
    for b in range(0, T, 2):                     # deletion indicator active for half of the targets
        DEL[rng.choice(n, 3 + b % 5, replace=False), b] = 1.0
    Y = rng.normal(size=(n, T)) + 0.2 * X[:, :T % R + 1].sum(1, keepdims=True) + 0.5 * XB + 2 * DEL
    return D, X, XB, DEL, Y


@pytest.mark.parametrize("seed", [0, 1])
def test_closed_form_vs_ols_toy(seed):
    D, X, XB, DEL, Y = _toy(seed)
    res = Residualiser(D)
    assert res.rank == D.shape[1] - 1                       # rank-deficient design detected
    r, t, df = pair_partial(res(X), res(Y), [], [res(XB), res(DEL)], df_base=res.df_base)
    n_checked = 0
    for a in range(X.shape[1]):
        for b in range(Y.shape[1]):
            cols = [D, XB[:, b]] + ([DEL[:, b]] if DEL[:, b].any() else []) + [X[:, a]]
            t_ols, df_ols = _ols_t(Y[:, b], cols)
            assert df[a, b] == df_ols                      # per-target df (k_B = 1 or 2)
            assert abs(t[a, b] - t_ols) <= TOL * max(1.0, abs(t_ols))
            n_checked += 1
    assert n_checked >= 200
    # chunking does not change anything
    _, t2, _ = chunked_pair_partial(res(X), res(Y), [], [res(XB), res(DEL)], df_base=res.df_base, chunk=5)
    np.testing.assert_allclose(t2, t, rtol=1e-12)


def test_closed_form_vs_ols_regulator_nuisance_and_coexpression():
    """k = 3 pair-specific regressors: hinge (gain slope as V_A) and co-expression (x_g as V_A)."""
    D, X, XB, DEL, Y = _toy(3, R=6, T=10)
    res = Residualiser(D)
    xl, xg = np.minimum(X, 0), np.maximum(X, 0)
    _, t, df = pair_partial(res(xl), res(Y), [res(xg)], [res(XB), res(DEL)], df_base=res.df_base)
    for a in range(X.shape[1]):
        for b in range(Y.shape[1]):
            cols = [D, xg[:, a], XB[:, b]] + ([DEL[:, b]] if DEL[:, b].any() else []) + [xl[:, a]]
            t_ols, df_ols = _ols_t(Y[:, b], cols)
            assert df[a, b] == df_ols
            assert abs(t[a, b] - t_ols) <= TOL * max(1.0, abs(t_ols))


def test_diag_partial_matches_ols():
    D, X, XB, DEL, Y = _toy(4, R=8, T=8)
    res = Residualiser(D)
    _, t, df = diag_partial(res(X), res(Y), df_base=res.df_base)
    for j in range(8):
        t_ols, df_ols = _ols_t(Y[:, j], [D, X[:, j]])
        assert df == df_ols and abs(t[j] - t_ols) <= TOL * max(1, abs(t_ols))


def test_guards():
    D, X, XB, DEL, Y = _toy(5, R=4, T=4, rank_deficient=False)
    X[:, 0] = D[:, 1] * 2 + 1e-9 * np.random.default_rng(0).normal(size=D.shape[0])   # explained by D
    res = Residualiser(D)
    U = res(X)
    _, t, _ = pair_partial(U, res(Y), [], [res(XB)], df_base=res.df_base,
                           u_raw=np.linalg.norm(X - X.mean(0), axis=0), y_raw=np.linalg.norm(Y - Y.mean(0), axis=0))
    assert np.isnan(t[0]).all() and np.isfinite(t[1:]).all()
    with pytest.raises(AssertionError):
        Residualiser(np.column_stack([D, D[:, 1] + 1e-9 * np.arange(D.shape[0])]))   # cond >= 1e6


def test_z_from_t_and_lambda():
    t = np.array([-60.0, -2.0, 0.0, 2.0, 60.0])
    z = z_from_t(t, 1000)
    assert np.all(np.isfinite(z)) and np.all(np.sign(z) == np.sign(t)) and 30 < abs(z[-1]) < 60 and abs(z[1]) < 2
    rng = np.random.default_rng(0)
    assert abs(lambda_gc(rng.standard_t(1300, 200_000), 1300) - 1) < 0.02


# ----------------------------------------------------------------------------- real data
real = pytest.mark.skipif(not (PROC / "panel_s2.tsv").exists(), reason="real panel absent")


@pytest.fixture(scope="module")
def real_prep():
    from depmaptr import io as dio
    from depmaptr.pipeline import Params, prepare
    ds = dio.load_real("S2", checksums=False)
    return prepare(ds, Params())


@real
@pytest.mark.parametrize("arm_handling", ["loco", "focal"])
def test_closed_form_vs_ols_real_pairs(real_prep, arm_handling):
    """>= 200 random real trans pairs (half with an active deletion indicator), production code path."""
    from depmaptr.covariates import DesignBuilder
    from depmaptr.pipeline import GroupCache, group_t
    prep = real_prep
    ds = prep.ds
    b = DesignBuilder(ds, prep.Z, prep.arms, arm_handling)
    prep.builder = b
    cache = GroupCache(prep, b)
    rng = np.random.default_rng(20261027)
    regs = np.asarray(ds.regulators)
    tch = ds.genes.loc[ds.targets, "chrom"].astype(str).to_numpy()
    del_active = np.flatnonzero(prep.DEL.any(0))
    checked = 0
    for c in ["1", "8", "17", "X"]:
        rk = regs[ds.genes.loc[regs, "chrom"].astype(str).to_numpy() == c]
        rk = rng.choice(rk, 6, replace=False)
        g = b.group_of(rk)[0]
        gv = cache.get(g)
        t, df, _ = group_t(prep, gv, rk)
        D, _ = b.design(g)
        cand_del = del_active[tch[del_active] != c]
        cand = np.flatnonzero(tch != c)
        tj = np.concatenate([rng.choice(cand_del, 5, replace=False), rng.choice(cand, 5, replace=False)])
        for i, a in enumerate(rk):
            xa = prep.x_all[:, prep.cn_idx[a]]
            for j in tj:
                cols = [D, prep.XT[:, j]] + ([prep.DEL[:, j]] if prep.DEL[:, j].any() else []) + [xa]
                t_ols, df_ols = _ols_t(prep.YT[:, j], cols)
                assert df[i, j] == df_ols, (a, ds.targets[j])
                assert abs(t[i, j] - t_ols) <= TOL * max(1.0, abs(t_ols)), (a, ds.targets[j], t[i, j], t_ols)
                checked += 1
    assert checked >= 200


@real
def test_saved_m1_scores_match_engine(real_prep):
    """If an M1 score file exists, its float32 t equals the engine t for a few rows (assembly check)."""
    from depmaptr.covariates import DesignBuilder
    from depmaptr.pipeline import GroupCache, group_t
    f = ROOT / "data" / "results" / "scores" / "M1" / "scores.npz"
    if not f.exists():
        pytest.skip("no M1 score file")
    z = np.load(f)
    prep = real_prep
    prep.builder = DesignBuilder(prep.ds, prep.Z, prep.arms, "loco")
    regs = list(z["regulators"])
    assert regs == list(prep.ds.regulators) and list(z["targets"]) == list(prep.ds.targets)
    rk = [regs[0], regs[5000], regs[-1]]
    for a in rk:
        g = prep.builder.group_of([a])[0]
        t, _, _ = group_t(prep, GroupCache(prep).get(g), [a])
        np.testing.assert_allclose(z["t"][regs.index(a)], t[0].astype(np.float32), rtol=1e-6, atol=1e-6)
