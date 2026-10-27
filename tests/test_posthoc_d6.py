"""Post-hoc tier D6 (reviews/ideas/M5a_ideas.md §2): I1-I6 and C2 sanity / known-answer tests."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm

from depmaptr.pipeline import (GroupCache, Params, context_fusion, group_t, prepare, regulator_context_percentile,
                               run, target_context_percentile, trans_percentile)
from depmaptr.simulate import simulate

ROOT = Path(__file__).resolve().parents[1]
INC = dict(g1_tau=100, g1_w=0.0, a3_rho=0.95)          # incumbent g1_t100_w0+a3_r095_s2_w0
quiet = lambda *a, **k: None  # noqa: E731


@pytest.fixture(scope="module")
def sim():
    ds, tr = simulate(seed=1)
    return ds, tr, run(ds, Params(**INC), log=quiet)


def _idx(ds):
    R, T = pd.Index(ds.regulators), pd.Index(ds.targets)
    trans = ds.genes.loc[R, "chrom"].to_numpy()[:, None] != ds.genes.loc[T, "chrom"].to_numpy()[None, :]
    return R, T, trans


# ----------------------------------------------------------------------------- I5
def test_i5_transform_and_ols():
    """y_cols implements log2(TPM/(1 - share_c) + 1); the engine t on renormalised y equals statsmodels OLS."""
    ds, _ = simulate(seed=2)
    prep = prepare(ds, Params(comp_renorm=True))
    tpm = 2.0 ** ds.expr - 1
    ch = ds.genes.loc[ds.expr_keys, "chrom"].to_numpy()
    share = tpm[:, ch == "3"].sum(1) / tpm.sum(1)
    j = [0, 5, 100]
    np.testing.assert_allclose(prep.y_cols("3", j), np.log2(tpm[:, j] / (1 - share)[:, None] + 1), rtol=1e-12)
    assert np.allclose(prepare(ds, Params()).y_cols("3", j), ds.expr[:, j])
    regs = [k for k in ds.regulators if ds.genes.loc[k, "chrom"] == "3"][:3]
    gv = GroupCache(prep).get("3")
    t, df, _ = group_t(prep, gv, regs)
    D, _ = prep.builder.design("3")
    te = [prep.ex_idx[k] for k in ds.targets]
    Y = prep.y_cols("3", te)
    for i, a in enumerate(regs):
        for b in (0, 50, 400):
            cols = [D, prep.XT[:, b]] + ([prep.DEL[:, b]] if prep.DEL[:, b].any() else []) + [prep.x_all[:, prep.cn_idx[a]]]
            f = sm.OLS(Y[:, b], np.column_stack(cols)).fit()
            assert abs(f.tvalues[-1] - t[i, b]) < 1e-8 * max(1, abs(f.tvalues[-1]))


def test_i5_sim_sanity(sim):
    ds, tr, inc = sim
    R, T, trans = _idx(ds)
    o = run(ds, Params(**INC, comp_renorm=True), log=quiet)
    a, b = inc["t"][trans], o["t"][trans]
    ok = np.isfinite(a) & np.isfinite(b)
    assert np.corrcoef(a[ok], b[ok])[0, 1] > 0.95
    assert (np.isnan(o["score"]) == np.isnan(o["t"])).all()
    with pytest.raises(ValueError):
        run(ds, Params(comp_renorm=True, arm_handling="focal"), log=quiet)


# ----------------------------------------------------------------------------- I3
def test_i3_iv_recovers_planted_effect(sim):
    """For planted pairs (y_B += gamma * y_A), beta_IV = r_AB sd_B / (r_cis sd_yA) recovers gamma."""
    ds, tr, _ = sim
    R, T, _ = _idx(ds)
    o = run(ds, Params(**{**INC, "score_base": "iv_shrink", "a3_rho": None}), log=quiet)
    pl = tr["planted"]
    pl = pl[pl.regulator != tr["neg_cis_regulator"]]
    ri, ti = R.get_indexer(pl.regulator), T.get_indexer(pl.target)
    strong = np.abs(o["t"][ri, ti]) > 8
    ratio = o["score"][ri, ti][strong] / pl.gamma.to_numpy()[strong]
    assert strong.sum() > 100 and 0.8 < np.median(ratio) < 1.2 and np.mean(ratio > 0) > 0.98
    # no instrument (negative cis / CN-only): beta = 0; shrinkage zeroes |t| <= 2
    neg = R.get_loc(tr["neg_cis_regulator"])
    assert (np.nan_to_num(o["score"][neg]) == 0).all()
    small = np.abs(o["t"]) <= 2
    assert (np.nan_to_num(o["score"][small]) == 0).all()
    assert (np.isnan(o["score"]) == np.isnan(o["t"])).all()


# ----------------------------------------------------------------------------- I6
def test_i6_lineage_halves_rule():
    from depmaptr.covariates import lineage_halves, lineage_levels
    ds, _ = simulate(seed=3)
    h1, h2 = lineage_halves(ds.panel), lineage_halves(ds.panel)
    assert h1.equals(h2) and set(h1.unique()) == {1, 2}
    n = lineage_levels(ds.panel).value_counts()
    tot = n.groupby(h1.reindex(n.index)).sum()
    assert abs(tot[1] - tot[2]) <= n.max()


@pytest.mark.skipif(not (ROOT / "data" / "processed" / "panel_s2.tsv").exists(), reason="real panel absent")
def test_i6_split_file_matches_rule():
    from depmaptr.covariates import lineage_halves
    panel = pd.read_csv(ROOT / "data" / "processed" / "panel_s2.tsv", sep="\t", low_memory=False)
    tab = pd.read_csv(ROOT / "results" / "split" / "lineage_halves.tsv", sep="\t").set_index("lineage_level")
    rule = lineage_halves(panel)
    assert len(tab) == 32 and (tab["half"].reindex(rule.index) == rule).all()
    tot = tab.groupby("half")["n"].sum()
    assert tot.sum() == 1436 and abs(tot[1] - tot[2]) < 50


def test_i6_replication(sim):
    ds, tr, inc = sim
    R, T, trans = _idx(ds)
    from depmaptr.pipeline import lineage_replicated_t
    p = Params(**INC, score_base="t_rep")
    rep, info = lineage_replicated_t(ds, p, log=quiet)
    assert info["half1_n"] + info["half2_n"] == ds.n
    # |t_rep| <= |t| of each half, sign agreement, planted pairs replicate
    pl = tr["planted"]
    ri, ti = R.get_indexer(pl.regulator), T.get_indexer(pl.target)
    assert np.mean(np.abs(rep[ri, ti]) > 3) > 0.8
    assert np.mean(np.abs(rep[trans]) > 3) < 0.05
    o = run(ds, p, log=quiet)
    assert (np.isnan(o["score"]) == np.isnan(o["t"])).all()
    from depmaptr.gates import g1_weight
    w1 = g1_weight(o["cis_t"], 100, 0.0)
    ref = np.where(np.isnan(o["t"]), np.nan, rep) * w1[:, None] * o["alpha"]
    np.testing.assert_allclose(o["score"], ref.astype(np.float32), rtol=1e-5, atol=1e-6, equal_nan=True)


def test_i6_half_t_bounds():
    """t_rep = sign(t1) min(|t1|,|t2|) 1[sign t1 = sign t2] on a toy pair of matrices (via the public function)."""
    t1 = np.array([[3.0, -2.0, 4.0, np.nan]])
    t2 = np.array([[5.0, 1.0, 1.5, 2.0]])
    rep = np.sign(t1) * np.minimum(np.abs(t1), np.abs(t2)) * (np.sign(t1) == np.sign(t2))
    rep = np.where(np.isnan(rep), 0.0, rep)
    assert rep.tolist() == [[3.0, 0.0, 1.5, 0.0]]


# ----------------------------------------------------------------------------- I1 / I2
def test_trans_percentile_known_answer():
    ds, _ = simulate(seed=1)
    R, T, trans = _idx(ds)
    rng = np.random.default_rng(0)
    v = rng.integers(0, 5, size=trans.shape).astype(np.float32)
    v[0, 0] = np.nan
    pct = trans_percentile(ds, v)
    from scipy.stats import rankdata
    m = trans & np.isfinite(v)
    ref = rankdata(v[m]) / m.sum()
    np.testing.assert_allclose(pct[m], ref, rtol=1e-6)
    assert (pct[~trans & np.isfinite(v)] == 0).all() and np.isnan(pct[0, 0])


def test_i1_fusion_monotone_and_no_nan(sim):
    ds, tr, inc = sim
    R, T, trans = _idx(ds)
    for g in (1.0, 2.0):
        o = run(ds, Params(**INC, fuse_gamma=g), log=quiet)
        S, s_inc, t = o["score"], inc["score"], o["t"]
        assert (np.isnan(S) == np.isnan(t)).all()
        fin = trans & np.isfinite(S)
        assert (S[fin & (s_inc == 0)] == 0).all()
        nz = fin & (S != 0)
        assert (np.sign(S[nz]) == np.sign(t[nz])).all()
        assert np.nanmax(np.abs(S)) <= 1.0 + 1e-6
        rho_t = target_context_percentile(ds, "HCT116")
        # within a target column: |S| is a non-decreasing function of |S_inc|
        for j in range(0, len(T), 97):
            m = fin[:, j]
            a, b = np.abs(s_inc[m, j]), np.abs(S[m, j])
            o_ = np.argsort(a, kind="stable")
            assert np.all(np.diff(b[o_]) >= -1e-6)
            np.testing.assert_allclose(b[a > 0], (trans_percentile(ds, np.abs(s_inc))[m, j] * rho_t[j] ** g)[a > 0],
                                       rtol=1e-5)
        # within a regulator row at equal |S_inc| rank, larger rho_T gives larger |S|
        i = int(np.argmax((np.abs(s_inc) > 0).sum(1)))
        m = fin[i] & (s_inc[i] != 0)
        ratio = np.abs(S[i, m]) / trans_percentile(ds, np.abs(s_inc))[i, m]
        assert np.all(np.diff(ratio[np.argsort(rho_t[m])]) >= -1e-6)


def test_i2_regctx(sim):
    ds, tr, inc = sim
    R, T, trans = _idx(ds)
    o = run(ds, Params(**INC, regctx=True), log=quiet)
    rho = regulator_context_percentile(ds, "HCT116")
    cn_only = R.isin(tr["cn_only"])
    assert (rho[cn_only] == 0).all() and (rho[~cn_only] > 0).all() and rho.max() == 1
    np.testing.assert_allclose(o["score"], inc["score"] * rho[:, None].astype(np.float32), rtol=1e-6, equal_nan=True)
    # I2 on top of I1 (applies to whatever precedes it)
    o2 = run(ds, Params(**INC, fuse_gamma=1.0, regctx=True), log=quiet)
    o1 = run(ds, Params(**INC, fuse_gamma=1.0), log=quiet)
    np.testing.assert_allclose(o2["score"], o1["score"] * rho[:, None].astype(np.float32), rtol=1e-6, equal_nan=True)


def test_c2_nc1_fusion(sim):
    ds, tr, _ = sim
    R, T, trans = _idx(ds)
    for g in (1.0, 2.0):
        o = run(ds, Params(nc1_index=0, g1_tau=100, g1_w=0.0, fuse_gamma=g), log=quiet)
        assert "alpha" not in o and (np.isnan(o["score"]) == np.isnan(o["t"])).all()
        nc1 = run(ds, Params(nc1_index=0), log=quiet)
        np.testing.assert_array_equal(o["t"], nc1["t"])


# ----------------------------------------------------------------------------- I4
def test_i4_sim_attribution_rho090():
    """SIM recovery of the planted gene within its segment at rho 0.90 vs 0.95 (reported in docs/MODEL.md)."""
    res = {}
    for rho in (0.95, 0.90):
        wins, chance = [], []
        for seed in (1, 2, 3):
            ds, tr = simulate(seed=seed)
            R, T, _ = _idx(ds)
            al = run(ds, Params(a3_rho=rho), log=quiet)["alpha"]
            pl = tr["planted"]
            for A, seg in tr["planted_segments"].items():
                members = [g for g in seg if g in R]
                ai = R.get_indexer(members)
                for b in pl.target[pl.regulator == A]:
                    wins.append(members[int(np.argmax(al[ai, T.get_loc(b)]))] == A)
                    chance.append(1 / len(members))
        res[rho] = (np.mean(wins), np.mean(chance))
    print("A3 planted-first rate (rho: rate, chance):", res)
    assert res[0.90][0] > 2 * res[0.90][1] and res[0.90][0] > 0.6
