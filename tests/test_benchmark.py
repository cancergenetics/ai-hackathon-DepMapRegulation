"""Known-answer tests for the benchmark harness (src/depmaptr/benchmark.py) on toy inputs.

Every metric is checked against a hand calculation or a dense brute-force reference."""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from depmaptr import benchmark as bm

NaN = np.nan


# ----------------------------------------------------------------------------- toy builders
def toy_tables(nr=40, nt=150, n_chrom=4, arms_per_chrom=2, rate=None, seed=0, n_bins_per_chrom=5):
    """Random toy universe: regulators/targets on `n_chrom` chromosomes, labels Bernoulli(rate_i)."""
    rng = np.random.default_rng(seed)
    chroms = [str(c + 1) for c in range(n_chrom)]
    rc = rng.choice(chroms, nr)
    arm = [c + ("p" if rng.random() < 0.5 or arms_per_chrom == 1 else "q") for c in rc]
    regs = pd.DataFrame({"regulator": [f"R{i}" for i in range(nr)], "chrom": rc, "arm": arm})
    regs["own_class"] = rng.choice(["i", "ii", "iii"], nr)
    regs["cn_tercile"] = rng.integers(1, 4, nr)
    regs["degree_decile"] = rng.integers(1, 11, nr)
    regs["fold"] = rng.integers(1, 3, nr)
    tc = rng.choice(chroms, nt)
    tgts = pd.DataFrame({"target": [f"T{j}" for j in range(nt)], "chrom": tc})
    tgts["bin"] = tgts["chrom"] + ":" + pd.Series(rng.integers(0, n_bins_per_chrom, nt)).astype(str)
    tgts["expr_decile"] = bm.quantile_bins(rng.random(nt), 10)
    if rate is None:
        rate = rng.uniform(0.0, 0.15, nr)
    L = rng.random((nr, nt)) < np.asarray(rate)[:, None]
    L &= rc[:, None] != tc[None, :]
    return regs, tgts, L


def pos_from_dense(regs, tgts, L, lfc=None, padj=None):
    r, t = np.nonzero(L)
    n = len(r)
    lfc = np.ones(n, np.int8) if lfc is None else lfc
    padj = np.full(n, 0.001) if padj is None else padj
    return pd.DataFrame({"regulator": regs["regulator"].values[r], "target": tgts["target"].values[t],
                         "p_adj": padj, "lfc_sign": lfc, "ueffect": 0.9,
                         "pos_primary": padj < 0.05, "pos_ueffect": padj < 0.05, "pos_p01": padj < 0.01})


def dense_valid(regs, tgts):
    return regs["chrom"].values[:, None] != tgts["chrom"].values[None, :]


def dense_sel(ev, frac):
    m = np.zeros(ev.nr * ev.nt, bool)
    m[ev.topk(frac)] = True
    return m.reshape(ev.nr, ev.nt)


def dense_or(sel, L, V, wr=None, wt=None):
    wr = np.ones(sel.shape[0]) if wr is None else wr
    wt = np.ones(sel.shape[1]) if wt is None else wt
    W = wr[:, None] * wt[None, :] * V
    a = (W * (sel & L)).sum()
    s = (W * sel).sum()
    p = (W * L).sum()
    n = W.sum()
    return bm.table_from_counts(a, s, p, n)


# ----------------------------------------------------------------------------- hand-computed OR
def hand_toy():
    regs = pd.DataFrame({"regulator": ["R1", "R2"], "chrom": ["1", "2"], "arm": ["1p", "2q"],
                         "own_class": ["i", "iii"], "cn_tercile": [1, 2], "degree_decile": [1, 2], "fold": [1, 2]})
    tgts = pd.DataFrame({"target": ["T1", "T2", "T3", "T4", "T5"], "chrom": ["1", "2", "3", "3", "1"],
                         "bin": ["1:0", "2:0", "3:0", "3:1", "1:1"], "expr_decile": [1, 2, 3, 4, 5]})
    S = np.array([[100, 5, -4, 1, 50],          # R1: T1, T5 same chromosome -> masked
                  [3, 99, NaN, 0.5, 2]], dtype=np.float32)   # R2: T2 masked, T3 NaN
    pos = pd.DataFrame({"regulator": ["R1", "R2", "R2"], "target": ["T3", "T5", "T4"], "p_adj": [0.01, 0.02, 0.001],
                        "lfc_sign": [-1, 1, 1], "ueffect": [0.1, 0.5, 0.9]})
    pos["pos_primary"] = True
    pos["pos_ueffect"] = pos["ueffect"].abs() >= 0.2
    pos["pos_p01"] = pos["p_adj"] < 0.01
    return regs, tgts, pos, S


def test_hand_or_and_masking():
    regs, tgts, pos, S = hand_toy()
    ev = bm.Evaluation(regs, tgts, pos, S, np.array([1, 1], np.int8))
    # R1 (chr1): T2, T3, T4 (3 pairs); R2 (chr2): T1, T3, T4, T5 (4 pairs)
    assert list(ev.n_i) == [3, 4] and ev.N == 7
    sel = ev.topk(3 / 7)                       # K = 3: R1T2 (5), R1T3 (4), R2T1 (3); masked 100/50/99 never
    names = {(ev.regs.regulator[i // 5], ev.tgts.target[i % 5]) for i in sel}
    assert names == {("R1", "T2"), ("R1", "T3"), ("R2", "T1")}
    r = ev.pooled(sel, ev.pos_lin())
    assert (r["a"], r["b"], r["c"], r["d"]) == (1, 2, 2, 2)
    assert r["or"] == pytest.approx(0.5)
    f = bm.fisher(1, 2, 2, 2)
    assert f["fisher_p_greater"] == pytest.approx(stats.fisher_exact([[1, 2], [2, 2]], "greater")[1])
    # NaN -> minimum: with K = 6 every non-NaN universe pair is selected and the NaN pair is not
    sel6 = ev.topk(6 / 7)
    assert 1 * 5 + 2 not in set(sel6)           # (R2, T3)
    assert len(sel6) == 6
    # sensitivity flags
    assert ev.pooled(sel, ev.pos_lin("pos_p01"))["c"] == 1   # only R2T4 (0.001) is p<0.01, unselected


def test_same_chromosome_never_counted():
    regs, tgts, L = toy_tables(nr=30, nt=80, seed=3)
    V = dense_valid(regs, tgts)
    S = np.where(~V, 1e6, np.random.default_rng(1).normal(size=V.shape)).astype(np.float32)
    ev = bm.Evaluation(regs, tgts, pos_from_dense(regs, tgts, L), S, np.ones(30, np.int8))
    assert ev.N == V.sum()
    for frac in (0.01, 0.1, 0.5, 1.0):
        sel = ev.topk(frac)
        r, t = np.divmod(sel, ev.nt)
        assert (regs["chrom"].values[r] != tgts["chrom"].values[t]).all()
    # all universe pairs selected -> table has b = n - p etc. and nothing same-chrom
    r = ev.pooled(ev.topk(1.0), ev.pos_lin())
    assert r["a"] + r["b"] == V.sum() and r["c"] == 0 and r["d"] == 0
    # a positives file containing a same-chromosome pair is rejected
    bad = pos_from_dense(regs, tgts, ~V & (np.arange(80)[None, :] == np.argmax(~V, 1)[:, None]))
    with pytest.raises(ValueError, match="same-chromosome"):
        bm.Evaluation(regs, tgts, bad, S, np.ones(30, np.int8))


def test_build_positives_drops_cis_and_self():
    regs = pd.DataFrame({"regulator": ["A (1)", "B (2)"], "pie_id": ["A", "B"], "chrom": ["1", "2"]})
    tgts = pd.DataFrame({"target": ["A (1)", "C (3)", "D (4)"], "pie_id": ["A", "C", "D"], "chrom": ["1", "1", "3"]})
    de = pd.DataFrame({"target": ["A", "A", "A", "B", "B", "B"], "feature": ["A", "C", "D", "C", "D", "Z"],
                       "log2_fold_change": [-25, 1, -1, 0.5, 2, 1], "p_adj": [1e-9, 1e-3, 0.03, 0.06, 0.004, 1e-4],
                       "Ueffect": [-0.9, 0.7, -0.05, 0.6, 0.3, 0.5]})
    p = bm.build_positives(de, regs, tgts)
    assert list(zip(p.regulator, p.target)) == [("A (1)", "D (4)"), ("B (2)", "D (4)")]
    assert list(p.pos_ueffect) == [False, True] and list(p.pos_p01) == [False, True]
    assert list(p.lfc_sign) == [-1, 1]


def test_tie_breaks_seeded_and_deterministic():
    regs, tgts, L = toy_tables(nr=20, nt=50, seed=5)
    S = np.ones((20, 50), np.float32)
    pos = pos_from_dense(regs, tgts, L)
    a = bm.Evaluation(regs, tgts, pos, S, np.ones(20, np.int8)).topk(0.1)
    b = bm.Evaluation(regs, tgts, pos, S, np.ones(20, np.int8)).topk(0.1)
    c = bm.Evaluation(regs, tgts, pos, S, np.ones(20, np.int8), seed=7).topk(0.1)
    assert np.array_equal(a, b) and not np.array_equal(a, c)
    # row order does not change which (regulator, target) pairs win
    perm = np.random.default_rng(0).permutation(20)
    ev2 = bm.Evaluation(regs.iloc[perm], tgts, pos, S, np.ones(20, np.int8))
    names = lambda ev, s: {(ev.regs.regulator[i // ev.nt], ev.tgts.target[i % ev.nt]) for i in s}
    ev1 = bm.Evaluation(regs, tgts, pos, S, np.ones(20, np.int8))
    assert names(ev1, ev1.topk(0.1)) == names(ev2, ev2.topk(0.1))


# ----------------------------------------------------------------------------- permutation
def test_permutation_fast_counts_match_bruteforce():
    regs, tgts, L = toy_tables(nr=25, nt=60, seed=11)
    V = dense_valid(regs, tgts)
    S = np.random.default_rng(2).normal(size=V.shape).astype(np.float32) + 2 * L
    ev = bm.Evaluation(regs, tgts, pos_from_dense(regs, tgts, L), S, np.ones(25, np.int8))
    sel = ev.topk(0.1)
    Sd = dense_sel(ev, 0.1)
    rng = np.random.default_rng(4)
    perms = [rng.permutation(25) for _ in range(20)] + [np.arange(25)]
    res = ev.permutation(sel, ev.pos_lin(), np.zeros(25), perms=perms, return_null=True)
    for j, (tn, to) in zip(perms, res["tables"]):
        P = V & V[j]                                   # trans for both A and pi(A)
        an = (Sd & L[j] & P).sum(); s = (Sd & P).sum(); pn = (L[j] & P).sum(); n = P.sum()
        assert tuple(map(int, tn)) == tuple(map(int, bm.table_from_counts(an, s, pn, n)))
        ao = (Sd & L & P).sum(); po = (L & P).sum()
        assert tuple(map(int, to)) == tuple(map(int, bm.table_from_counts(ao, s, po, n)))
    # identity permutation: null == observed == full OR
    assert res["null"][-1] == pytest.approx(res["obs"][-1]) == pytest.approx(res["or"])


def test_permutation_respects_strata():
    regs, tgts, L = toy_tables(nr=30, nt=60, seed=12)
    S = np.random.default_rng(0).normal(size=L.shape).astype(np.float32)
    ev = bm.Evaluation(regs, tgts, pos_from_dense(regs, tgts, L), S, np.ones(30, np.int8))
    strata = np.arange(30) % 3
    groups = ev._strata_index(strata)
    rng = np.random.default_rng(0)
    for _ in range(5):
        j = np.empty(30, int)
        for g in groups:
            j[g] = g[rng.permutation(len(g))]
        assert (strata[j] == strata).all()
    # singleton strata -> every permutation is the identity -> p = 1
    r = ev.permutation(ev.topk(0.1), ev.pos_lin(), np.arange(30), n_perm=19)
    assert r["perm_p"] == 1.0


def test_permutation_strong_signal_minimum_p():
    regs, tgts, L = toy_tables(nr=40, nt=100, seed=13, rate=np.full(40, 0.08))
    S = (np.random.default_rng(1).normal(size=L.shape) + 4 * L).astype(np.float32)
    ev = bm.Evaluation(regs, tgts, pos_from_dense(regs, tgts, L), S, np.ones(40, np.int8))
    r = ev.permutation(ev.topk(0.05), ev.pos_lin(), np.zeros(40), n_perm=99, rng=np.random.default_rng(0))
    assert r["perm_p"] == pytest.approx(0.01) and r["or"] > 10 and r["or_over_median_null"] > 5


@pytest.mark.parametrize("which", ["regulator", "degree"])
def test_permutation_calibrated_under_random_scores(which):
    """B0-style null: random scores -> permutation p approximately uniform over 200 seeds.

    The toy has 23 chromosomes like the real universe. With very few chromosomes (e.g. 4) and tiny counts
    the test is conservative (p < 0.05 in ~1 % of null runs), because each replicate's pair set then drops
    a large share of the pairs; see docs/BENCHMARK.md."""
    ps = []
    for seed in range(200):
        regs, tgts, L = toy_tables(nr=60, nt=400, n_chrom=23, seed=1000 + seed)
        S = np.random.default_rng(1000 + seed).normal(size=L.shape).astype(np.float32)
        ev = bm.Evaluation(regs, tgts, pos_from_dense(regs, tgts, L), S, np.ones(60, np.int8))
        if which == "regulator":       # own-gene class x CN tercile
            strata = (regs.own_class.astype(str) + "|" + regs.cn_tercile.astype(str)).values
        else:                          # own-gene class x benchmark-degree decile
            strata = (regs.own_class.astype(str) + "|" + pd.Series(bm.quantile_bins(L.sum(1), 10)).astype(str)).values
        r = ev.permutation(ev.topk(0.05), ev.pos_lin(), strata, n_perm=99, rng=np.random.default_rng(seed))
        ps.append(r["perm_p"])
    assert_uniformish(ps)


def assert_uniformish(ps):
    """~uniform p over 200 null seeds: omnibus KS plus tail rates (binomial SE at 0.05 is ~0.015)."""
    ps = np.asarray(ps)
    assert stats.kstest(ps, "uniform").pvalue > 1e-3
    assert 0.01 <= np.mean(ps <= 0.05) <= 0.10
    assert 0.04 <= np.mean(ps <= 0.10) <= 0.17
    assert abs(np.mean(ps) - 0.5) < 0.06


# ----------------------------------------------------------------------------- jackknife / bootstrap
def test_jackknife_matches_bruteforce_and_t_critical():
    regs, tgts, L = toy_tables(nr=40, nt=100, seed=21, n_chrom=4, arms_per_chrom=2)
    V = dense_valid(regs, tgts)
    S = (np.random.default_rng(3).normal(size=L.shape) + 1.5 * L).astype(np.float32)
    ev = bm.Evaluation(regs, tgts, pos_from_dense(regs, tgts, L), S, np.ones(40, np.int8))
    sel = ev.topk(0.05)
    Sd = dense_sel(ev, 0.05)
    jk = ev.jackknife([sel], [ev.pos_lin()])
    arms = sorted(regs.arm.unique(), key=bm.arm_sort_key)
    G = len(arms)
    th = []
    for a in arms:
        keep = (regs.arm != a).values
        th.append(bm.log_or(*dense_or(Sd[keep], L[keep], V[keep])))
    th = np.array(th)
    full = bm.log_or(*dense_or(Sd, L, V))
    se = math.sqrt((G - 1) / G * ((th - th.mean()) ** 2).sum())
    assert jk["G"] == G and jk["estimate"] == pytest.approx(full)
    assert jk["se"] == pytest.approx(se)
    assert jk["t_crit"] == pytest.approx(stats.t.ppf(0.975, G - 1))
    assert jk["ci95"][0] == pytest.approx(full - stats.t.ppf(0.975, G - 1) * se)
    # hand value of the critical t for G = 4 arms
    assert stats.t.ppf(0.975, 3) == pytest.approx(3.182446, abs=1e-6)
    # TRAIN folds: OR over each fold's regulators, top-K fixed at the full-split selection
    for f in (1, 2):
        k = (regs.fold == f).values
        assert ev.fold_log_or(sel, ev.pos_lin())[f"fold{f}"] == pytest.approx(
            bm.log_or(*dense_or(Sd[k], L[k], V[k])))
    # paired jackknife of identical scores: zero difference, zero SE
    pj = ev.jackknife([sel, sel], [ev.pos_lin(), ev.pos_lin()])
    assert pj["estimate"] == 0 and pj["se"] == 0


def test_bootstrap_replicate_matches_bruteforce():
    regs, tgts, L = toy_tables(nr=40, nt=100, seed=22)
    V = dense_valid(regs, tgts)
    S = (np.random.default_rng(3).normal(size=L.shape) + L).astype(np.float32)
    ev = bm.Evaluation(regs, tgts, pos_from_dense(regs, tgts, L), S, np.ones(40, np.int8))
    sel = ev.topk(0.05)
    Sd = dense_sel(ev, 0.05)
    ones = (np.ones((1, ev.G)), np.ones((1, ev.H)))
    assert ev.bootstrap_log_or(sel, ev.pos_lin(), ones)[0] == pytest.approx(bm.log_or(*dense_or(Sd, L, V)))
    Wa, Wb = ev.boot_weights(5, np.random.default_rng(9))
    rep = ev.bootstrap_log_or(sel, ev.pos_lin(), (Wa, Wb))
    for b in range(5):
        wr = Wa[b][ev.ga]
        wt = Wb[b][ev.tb]
        assert rep[b] == pytest.approx(bm.log_or(*dense_or(Sd, L, V, wr, wt)))
    # CMH bootstrap: unit weights reproduce the point estimate; a random replicate matches brute force
    cm = ev.cmh(sel, ev.pos_lin())["cmh_or"]
    assert ev.bootstrap_log_cmh(sel, ev.pos_lin(), ones)[0] == pytest.approx(math.log(cm))
    repc = ev.bootstrap_log_cmh(sel, ev.pos_lin(), (Wa, Wb), chunk=2)
    for b in range(5):
        wr, wt = Wa[b][ev.ga], Wb[b][ev.tb]
        num = den = 0.0
        for i in range(ev.nr):
            for e in range(1, 11):
                m = V[i] & (tgts.expr_decile.values == e)
                a, bb, c, d = dense_or(Sd[i:i + 1, m], L[i:i + 1, m], V[i:i + 1, m], None, wt[m])
                n = a + bb + c + d
                if n > 0:
                    num += wr[i] * a * d / n
                    den += wr[i] * bb * c / n
        assert repc[b] == pytest.approx(math.log(num / den))


# ----------------------------------------------------------------------------- CMH
def test_mantel_haenszel_by_hand_and_statsmodels():
    from statsmodels.stats.contingency_tables import StratifiedTable
    # two strata: (a,b,c,d) = (10,20,5,40) and (3,7,6,30)
    a, b, c, d = np.array([10, 3]), np.array([20, 7]), np.array([5, 6]), np.array([40, 30])
    n = a + b + c + d
    hand = (10 * 40 / 75 + 3 * 30 / 46) / (20 * 5 / 75 + 7 * 6 / 46)
    r = bm.mantel_haenszel(a, b, c, d)
    assert r["cmh_or"] == pytest.approx(hand)
    st = StratifiedTable([np.array([[a[k], b[k]], [c[k], d[k]]]) for k in range(2)])
    assert r["cmh_or"] == pytest.approx(st.oddsratio_pooled)
    assert r["cmh_log_se"] == pytest.approx(st.logodds_pooled_se)
    assert r["cmh_p"] == pytest.approx(st.test_null_odds(correction=False).pvalue)


def test_cmh_tables_match_bruteforce():
    from statsmodels.stats.contingency_tables import StratifiedTable
    regs, tgts, L = toy_tables(nr=20, nt=80, seed=31)
    V = dense_valid(regs, tgts)
    S = (np.random.default_rng(3).normal(size=L.shape) + L).astype(np.float32)
    ev = bm.Evaluation(regs, tgts, pos_from_dense(regs, tgts, L), S, np.ones(20, np.int8))
    sel = ev.topk(0.1)
    Sd = dense_sel(ev, 0.1)
    tabs = []
    for i in range(20):
        for e in range(1, 11):
            m = V[i] & (tgts.expr_decile.values == e)
            tabs.append(dense_or(Sd[i:i + 1, m], L[i:i + 1, m], V[i:i + 1, m]))
    tabs = np.array(tabs, float)
    tabs = tabs[tabs.sum(1) > 0]
    st = StratifiedTable([t.reshape(2, 2) for t in tabs])
    assert ev.cmh(sel, ev.pos_lin())["cmh_or"] == pytest.approx(st.oddsratio_pooled)
    # own-gene class x decile strata
    tabs = []
    for cls in sorted(regs.own_class.unique()):
        k = (regs.own_class == cls).values
        for e in range(1, 11):
            m = V[k] & (tgts.expr_decile.values == e)[None, :]
            tabs.append(dense_or(Sd[k] & m, L[k] & m, m))
    st = StratifiedTable([np.array(t, float).reshape(2, 2) for t in tabs if sum(t) > 0])
    assert ev.cmh(sel, ev.pos_lin(), "own_class")["cmh_or"] == pytest.approx(st.oddsratio_pooled)


# ----------------------------------------------------------------------------- sign metrics
def test_sign_concordance_c0_and_negative_cis():
    regs = pd.DataFrame({"regulator": ["R1", "R2"], "chrom": ["1", "2"], "arm": ["1p", "2p"],
                         "own_class": ["i", "i"], "cn_tercile": [1, 1], "degree_decile": [1, 1], "fold": [1, 2]})
    tgts = pd.DataFrame({"target": [f"T{j}" for j in range(1, 7)], "chrom": ["3"] * 6, "bin": ["3:0"] * 6,
                         "expr_decile": [1] * 6})
    S = np.full((2, 6), 0.1, np.float32)
    S[0, 0], S[0, 1], S[0, 4] = 5, -4, 6         # R1 (s_cis = +1): T1, T2, T5
    S[1, 2], S[1, 3] = 3, -2                     # R2 (s_cis = -1): T3, T4
    s_cis = np.array([1, -1], np.int8)
    pos = pd.DataFrame({"regulator": ["R1", "R1", "R2", "R2", "R1", "R1"],
                        "target": ["T1", "T2", "T3", "T4", "T5", "T6"],
                        "p_adj": 0.01, "lfc_sign": [-1, -1, 1, 1, -1, 1], "ueffect": 0.9})
    pos["pos_primary"] = pos["pos_ueffect"] = pos["pos_p01"] = True
    ev = bm.Evaluation(regs, tgts, pos, S, s_cis)
    sel = ev.topk(5 / 12)
    # predicted knockdown sign = -sign(score) * s_cis
    assert list(ev.pred_kd_sign(sel)) == [-1, 1, -1, 1, -1]     # R1T1, R1T2, R1T5, R2T3, R2T4
    m = ev.sign_metrics(sel, n_perm=50, min_overlap=2, rng=np.random.default_rng(0))
    # overlaps: pred [-1,+1,+1,-1,-1] vs lfc [-1,-1,+1,+1,-1] (R1T1, R1T2, R2T3, R2T4, R1T5)
    assert m["n_overlap"] == 5
    assert m["concordance"] == pytest.approx(3 / 5)
    assert m["q_pred_down"] == pytest.approx(3 / 5) and m["u_lfc_down"] == pytest.approx(3 / 5)
    assert m["c0"] == pytest.approx(0.6 * 0.6 + 0.4 * 0.4)
    assert m["sign_2x2"] == [2, 1, 1, 1] and m["sign_or"] == pytest.approx(2.0)
    assert m["regulator_weighted_concordance"] == pytest.approx((2 / 3 + 1 / 2) / 2)
    assert m["sign_perm_p"] == 1.0              # every regulator's LFC signs are constant -> null == obs
    # sign-concordant OR: concordant positives R1T1, R1T5, R2T3 (R1T6 pred down vs LFC up -> not)
    assert m["sign_concordant_table"] == [3, 2, 0, 7]
    assert m["sign_concordant_or"] == pytest.approx(3.5 * 7.5 / (2.5 * 0.5))


def test_sign_permutation_calibrated():
    """Random predicted and benchmark signs -> within-regulator sign-permutation p ~ uniform.
    (With only a handful of overlaps per regulator the discrete test is conservative, but still valid.)"""
    ps = []
    for seed in range(200):
        rng = np.random.default_rng(seed)
        regs, tgts, L = toy_tables(nr=40, nt=300, n_chrom=23, seed=seed, rate=np.full(40, 0.3))
        S = rng.normal(size=L.shape).astype(np.float32)
        pos = pos_from_dense(regs, tgts, L, lfc=rng.choice([-1, 1], L.sum()).astype(np.int8))
        ev = bm.Evaluation(regs, tgts, pos, S, rng.choice([-1, 1], 40).astype(np.int8))
        ps.append(ev.sign_metrics(ev.topk(0.3), n_perm=99, rng=rng)["sign_perm_p"])
    assert_uniformish(ps)


# ----------------------------------------------------------------------------- other metrics
def test_auroc_top10_callset_hub():
    regs, tgts, L = toy_tables(nr=30, nt=100, seed=41, rate=np.full(30, 0.2))
    V = dense_valid(regs, tgts)
    S = (np.random.default_rng(3).normal(size=L.shape) + L).astype(np.float32)
    S[0, :] = np.nan                                            # NaN row -> ties at the minimum
    ev = bm.Evaluation(regs, tgts, pos_from_dense(regs, tgts, L), S, np.ones(30, np.int8))
    pl = ev.pos_lin()
    au = ev.auroc_per_regulator(pl, min_pos=10)
    from sklearn.metrics import roc_auc_score
    ref = [roc_auc_score(L[i, V[i]], np.nan_to_num(np.abs(S[i, V[i]]), nan=-1)) for i in range(30)
           if (L[i] & V[i]).sum() >= 10]
    assert au["n_regulators"] == len(ref) and au["mean_auroc"] == pytest.approx(np.mean(ref))
    sel = ev.topk(0.1)
    cnt = np.sort(np.bincount(sel // ev.nt, minlength=30))[::-1]
    assert ev.top10_share(sel) == pytest.approx(cnt[:10].sum() / len(sel))
    calls = np.abs(np.nan_to_num(S)) > 2
    cs = ev.callset(calls)
    assert tuple(cs["table"]) == tuple(int(x) for x in dense_or(calls & V, L, V))
    hub = ev.hub_excluded(0.1)
    assert hub["n_reg_dropped"] == 3 and hub["n_tgt_dropped"] == 10


def test_quantile_bins_ties_share_bins():
    v = [0, 0, 0, 0, 0, 0, 1, 2, 3, 4]
    b = bm.quantile_bins(v, 10)
    assert len(set(b[:6])) == 1 and list(b[6:]) == sorted(b[6:])
    assert list(bm.quantile_bins(np.arange(9), 3)) == [1, 1, 1, 2, 2, 2, 3, 3, 3]
    assert bm.own_gene_class([True, True, True, False], [0.01, 0.2, 0.01, np.nan],
                             [-3, -3, 2, np.nan]).tolist() == ["i", "ii", "ii", "iii"]


# ----------------------------------------------------------------------------- split
def test_split_triplets_folds_and_seed():
    arms = [f"{c}{a}" for c in list(range(1, 23)) + ["X"] for a in "pq"]
    rng = np.random.default_rng(0)
    ac = pd.DataFrame({"arm": [str(a) for a in arms], "n_hct116": rng.integers(0, 300, len(arms)),
                       "n_hek293t": rng.integers(0, 100, len(arms))})
    ac.loc[0, ["n_hct116", "n_hek293t"]] = 0                      # 1p: no regulators -> not a unit
    ac.loc[1, ["n_hct116", "n_hek293t"]] = [0, 5]                 # 1q: HEK only -> a unit
    s = bm.make_split(ac)
    assert "1p" not in set(s.arm) and "1q" in set(s.arm)
    assert list(s.n_hct116) == sorted(s.n_hct116, reverse=True)
    draws = np.random.default_rng(20261027)
    for t, g in s.groupby("triplet"):
        k = int(draws.integers(3))
        exp = ["TEST" if i == k else "TRAIN" for i in range(len(g))]
        assert list(g.split) == exp
    tr = s[s.split == "TRAIN"]
    assert list(tr.fold) == [1 + i % 2 for i in range(len(tr))]
    assert (s[s.split == "TEST"].fold == 0).all()
    assert s.equals(bm.make_split(ac))


# ----------------------------------------------------------------------------- universe builders
def test_universe_builders():
    pie_map = pd.DataFrame({
        "pie_id": ["A", "B", "C", "D", "E", "F", "G", "H"],
        "depmap_cn_column": ["A (1)", "B (2)", None, "D (4)", "E (5)", "F (6)", "G (7)", "H (8)"],
        "depmap_expr_column": ["A (1)", None, None, "D (4)", "E (5)", "F (6)", "G (7)", None]})
    flags = pd.DataFrame({
        "gene": ["A (1)", "B (2)", "D (4)", "E (5)", "F (6)", "G (7)", "H (8)"],
        "chrom": ["1", "2", "X", "Y", "3", "X", "4"], "arm": ["p", "q", "p", "p", "q", "p", "p"],
        "variable_cn": [True, True, True, True, False, True, True],
        "measured": [True, False, True, True, True, True, False],
        "is_par": [False, False, False, False, False, True, False], "is_y": [False, False, False, True, False, False, False],
        "is_mt": False})
    regs, cnt = bm.build_universe_regulators(pie_map, flags, ["A", "B", "C", "D", "E", "F", "G"])
    assert list(regs.regulator) == ["A (1)", "B (2)", "D (4)"]     # C unmapped, E Y, F not variable, G PAR
    assert list(regs.arm) == ["1p", "2q", "Xp"]
    assert cnt["unmapped"] == 1 and cnt["excluded_Y_PAR_MT_unplaced"] == 2 and cnt["not_variable_cn"] == 1
    pos = pd.DataFrame({"start": [1e6, 0, 12e6, 0, 0, 0, 0], "end": [1.2e6, 0, 12.2e6, 0, 0, 0, 0]},
                       index=flags.gene)
    tg, tc = bm.build_universe_targets(pie_map, flags, ["A", "B", "D", "E", "H"], pos)
    assert list(tg.target) == ["A (1)", "D (4)"] and list(tg.bin) == ["1:0", "X:2"]
    assert tc["not_measured"] == 2
    bad = pie_map.copy()
    bad.loc[1, "depmap_cn_column"] = "A (1)"
    with pytest.raises(ValueError, match="many-to-one"):
        bm.pie_key_map(bad)


# ----------------------------------------------------------------------------- end-to-end + TEST guard
def write_toy_bench(root, seed=0):
    """Toy project tree: bench tables (TRAIN + TEST), sealed-style TEST labels, and two score files."""
    regs, tgts, L = toy_tables(nr=48, nt=120, n_chrom=4, arms_per_chrom=2, seed=seed)
    bd = root / "data" / "processed" / "bench"
    bd.mkdir(parents=True)
    arm_test = sorted(regs.arm.unique())[:2]
    regs["split"] = np.where(regs.arm.isin(arm_test), "TEST", "TRAIN")
    regs["fold"] = np.where(regs.split == "TRAIN", regs["fold"], 0)
    lfc = np.random.default_rng(seed).choice([-1, 1], L.sum()).astype(np.int8)
    pos = pos_from_dense(regs, tgts, L, lfc=lfc)
    test_regs = set(regs.regulator[regs.split == "TEST"])
    pos[~pos.regulator.isin(test_regs)].to_parquet(bd / "positives_hct116_TRAIN.parquet")
    pos[pos.regulator.isin(test_regs)].to_parquet(bd / "positives_hct116_TEST.parquet")
    lab = regs.loc[regs.split == "TEST", ["regulator", "own_class", "degree_decile"]]
    lab.to_csv(bd / "regulator_labels_hct116_TEST.tsv", sep="\t", index=False)
    r2 = regs.copy()
    r2.loc[r2.split == "TEST", ["own_class", "degree_decile"]] = np.nan
    r2.to_csv(bd / "universe_regulators_hct116.tsv", sep="\t", index=False)
    tgts.to_csv(bd / "universe_targets_hct116.tsv", sep="\t", index=False)
    runs = {}
    for name, strength in (("good", 2.0), ("random", 0.0)):
        rd = root / "data" / "results" / "scores" / name
        rd.mkdir(parents=True)
        S = (np.random.default_rng(seed + 1).normal(size=L.shape) + strength * L).astype(np.float32)
        # score file holds more regulators/targets than the universe, in another order
        extra_r = np.array(["XR1"]); extra_t = np.array(["XT1", "XT2"])
        R = np.concatenate([regs.regulator.values[::-1], extra_r]).astype(str)
        T = np.concatenate([extra_t, tgts.target.values]).astype(str)
        full = np.full((len(R), len(T)), 0.0, np.float32)
        full[:48, 2:] = S[::-1]
        np.savez(rd / "scores.npz", regulators=R, targets=T, t=full, score=full,
                 s_cis=np.ones(len(R), np.int8), cis_t=np.ones(len(R), np.float32))
        (rd / "meta.json").write_text(json.dumps({"run_id": name, "variant": name, "params": {"s": strength}}))
        runs[name] = rd
    return regs, tgts, L, runs


def test_evaluate_end_to_end_and_log(tmp_path):
    regs, tgts, L, runs = write_toy_bench(tmp_path)
    out = bm.evaluate(runs["good"], "hct116", "TRAIN", n_perm=49, n_boot=50, root=tmp_path)
    m = out["metrics"]
    tr = (regs.split == "TRAIN").values
    assert m["n_regulators"] == tr.sum() and m["n_targets"] == len(tgts)
    assert m["n_pairs"] == dense_valid(regs[tr], tgts).sum()
    assert m["primary_or"] > 3 and m["perm_regulator"]["perm_p"] < 0.1
    assert m["n_regulators_missing_in_scores"] == 0
    # same metrics as a direct Evaluation on the aligned dense block
    ev = bm.Evaluation(regs[tr], tgts, pos_from_dense(regs[tr], tgts, L[tr]),
                       (np.random.default_rng(1).normal(size=L.shape) + 2 * L).astype(np.float32)[tr],
                       np.ones(tr.sum(), np.int8))
    assert ev.pooled(ev.topk(bm.PRIMARY_K), ev.pos_lin())["or"] == pytest.approx(m["primary_or"])
    log = pd.read_csv(tmp_path / "results" / "tuning_log.tsv", sep="\t")
    assert len(log) == 1 and log.run_id[0] == "good" and log.split[0] == "TRAIN"
    for c in ("params_hash", "git_hash", "data_checksums", "seed", "metrics", "wall_s"):
        assert c in log.columns and pd.notna(log[c][0])
    assert json.loads(log.metrics[0])["primary_or"] == pytest.approx(m["primary_or"])
    assert (tmp_path / "results" / "eval" / "good" / "hct116_TRAIN.json").exists()
    bm.evaluate(runs["random"], "hct116", "TRAIN", n_perm=19, n_boot=20, root=tmp_path, full=False)
    assert len(pd.read_csv(tmp_path / "results" / "tuning_log.tsv", sep="\t")) == 2
    c = bm.compare(runs["good"], runs["random"], "hct116", "TRAIN", n_boot=50, root=tmp_path)
    assert c["delta_log_or"] > 0 and c["boot_delta_log_or_lb90"] > 0
    c0 = bm.compare(runs["good"], runs["good"], "hct116", "TRAIN", n_boot=20, root=tmp_path)
    assert c0["delta_log_or"] == 0 and c0["boot_delta_log_or_lb90"] == 0 and not c0["adopt_rule_8_2"]


def test_test_split_is_refused(tmp_path):
    _, _, _, runs = write_toy_bench(tmp_path)
    with pytest.raises(bm.SealedSplitError):
        bm.evaluate(runs["good"], "hct116", "TEST", root=tmp_path)
    with pytest.raises(bm.SealedSplitError):
        bm.compare(runs["good"], runs["random"], "hct116", "TEST", root=tmp_path)
    with pytest.raises(bm.SealedSplitError):
        bm.load_bench("hct116", "TEST", tmp_path)
    with pytest.raises(bm.SealedSplitError):
        bm.read_positives(tmp_path / "data/processed/bench/positives_hct116_TEST.parquet")
    with pytest.raises(bm.SealedSplitError):                # pytest is not scripts/final_eval.py
        bm.final_eval_token(tmp_path)
    with pytest.raises(bm.SealedSplitError):                # a forged object is not a token
        bm.evaluate(runs["good"], "hct116", "TEST", root=tmp_path, _token=object())
    assert not (tmp_path / "results" / "tuning_log.tsv").exists()


def test_final_eval_token_path_allows_test(tmp_path, monkeypatch):
    regs, tgts, L, runs = write_toy_bench(tmp_path)
    monkeypatch.setattr("sys.argv", [str(tmp_path / "scripts" / "final_eval.py")])
    tok = bm.final_eval_token(tmp_path / "data" / "processed" / "bench")
    out = bm.evaluate(runs["good"], "hct116", "TEST", n_perm=19, n_boot=20, root=tmp_path, _token=tok)
    assert out["split"] == "TEST" and out["metrics"]["n_regulators"] == (regs.split == "TEST").sum()
    mock = bm.evaluate(runs["good"], "hct116", "TEST", n_perm=19, n_boot=20, root=tmp_path, _token=tok,
                       _mock=True)
    assert mock["split"] == "PSEUDO_TEST"
    assert mock["metrics"]["n_regulators"] == ((regs.split == "TRAIN") & (regs.fold == 2)).sum()
