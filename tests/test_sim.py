"""SIM known-truth tests (DESIGN §7): recovery, confounder removal, attribution, sign rule, calibration."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

from depmaptr.baselines import run_baseline
from depmaptr.calibrate import lambda_gc
from depmaptr.pipeline import Params, run
from depmaptr.simulate import simulate

SEEDS = [1, 2, 3]
quiet = lambda *a, **k: None  # noqa: E731


def _index(ds):
    R, T = pd.Index(ds.regulators), pd.Index(ds.targets)
    rch = ds.genes.loc[R, "chrom"].to_numpy()
    tch = ds.genes.loc[T, "chrom"].to_numpy()
    return R, T, rch[:, None] != tch[None, :]


@pytest.fixture(scope="module", params=SEEDS)
def sim(request):
    ds, truth = simulate(seed=request.param)
    return ds, truth, run(ds, Params(), log=quiet), run_baseline(ds, "B1", log=quiet)


def test_planted_recovery_auroc(sim):
    ds, tr, m1, b1 = sim
    R, T, trans = _index(ds)
    pos = np.zeros(trans.shape, bool)
    pos[R.get_indexer(tr["planted"].regulator), T.get_indexer(tr["planted"].target)] = True
    neighbours = set(sum(tr["planted_segments"].values(), [])) - set(tr["planted_regulators"])
    reg_ok = ~R.isin(list(neighbours)) & (ds.genes.loc[R, "arm_id"].to_numpy() != tr["driver_arm"])
    ev = trans & reg_ok[:, None]
    auc = roc_auc_score(pos[ev], np.nan_to_num(np.abs(m1["score"][ev]), nan=0))
    assert auc > 0.98
    # B0 is at chance
    b0 = run_baseline(ds, "B0", log=quiet, seed=7)
    assert abs(roc_auc_score(pos[ev], b0["score"][ev]) - 0.5) < 0.1


# ----------------------------------------------------------------------------- ablation helpers
def _ablate(monkeypatch, what: str):
    """Remove one adjustment from the model (for negative-control 'the test must fail' checks)."""
    import depmaptr.pipeline as pl
    from depmaptr.covariates import DesignBuilder
    if what == "arms":
        orig = DesignBuilder.design

        def design(self, g):
            D, names = orig(self, g)
            keep = [i for i, nm in enumerate(names) if not nm.startswith("arm[")]
            return D[:, keep], [names[i] for i in keep]
        monkeypatch.setattr(DesignBuilder, "design", design)
    elif what in ("lineage", "sex"):
        orig_z = pl.build_Z

        def build_Z(panel):
            Z, log = orig_z(panel)
            Z = Z[[c for c in Z.columns if not c.startswith(what + "[")]]
            return Z, {**log, "columns": list(Z.columns)}
        monkeypatch.setattr(pl, "build_Z", build_Z)
    elif what == "homdel":
        monkeypatch.setattr(pl, "HOMDEL_MIN_LINES", 10 ** 9)
    else:
        raise ValueError(what)


def _fp(t, rows, cols, thr=3.0):
    return float(np.mean(np.abs(t[np.ix_(rows, cols)]) > thr))     # nominal rate at |t| > 3 is 0.0027


def test_cross_chromosome_arm_confounding(sim, monkeypatch):
    """1q and 3p move together (both directions, equal magnitude) and 3p dosage drives `driver_targets`:
    B1 and M1-without-arm-block report 1q -> driver-target false positives; M1 does not."""
    ds, tr, m1, b1 = sim
    R, T, _ = _index(ds)
    rows = ds.genes.loc[R, "arm_id"].to_numpy() == tr["coselect"][0]
    cols = T.isin(tr["driver_targets"])
    assert _fp(m1["t"], rows, cols) < 0.02
    assert _fp(b1["t"], rows, cols) > 0.1
    _ablate(monkeypatch, "arms")
    assert _fp(run(ds, Params(), log=quiet)["t"], rows, cols) > 0.1          # the adjustment is what works
    # the true LOCO effect of the driver arm is found by M1
    md = ds.genes.loc[R, "arm_id"].to_numpy() == tr["driver_arm"]
    assert _fp(m1["t"], md, cols) > 0.5


def test_lineage_confounding_removed(sim, monkeypatch):
    ds, tr, m1, b1 = sim
    R, T, _ = _index(ds)
    rows = ds.genes.loc[R, "arm_id"].to_numpy() == tr["lineage_arm"]
    cols = T.isin(tr["lineage_programme"])
    assert _fp(b1["t"], rows, cols) > 0.3 and _fp(m1["t"], rows, cols) < 0.02
    _ablate(monkeypatch, "lineage")
    assert _fp(run(ds, Params(), log=quiet)["t"], rows, cols) > 0.3


def test_sex_confounding_removed(sim, monkeypatch):
    """Raw X CN tracks sex and an autosomal programme is sex-dependent: X regulators must not hit it."""
    ds, tr, m1, b1 = sim
    R, T, _ = _index(ds)
    rows = (ds.genes.loc[R, "chrom"].to_numpy() == "X") & ~R.isin([tr["x_planted_regulator"]])
    cols = T.isin(tr["sex_programme"])
    assert _fp(m1["t"], rows, cols) < 0.02
    _ablate(monkeypatch, "sex")
    assert _fp(run(ds, Params(), log=quiet)["t"], rows, cols) > 0.3
    # the planted X regulator's targets are recovered by M1
    pl = tr["planted"]
    xr = tr["x_planted_regulator"]
    tt = m1["t"][R.get_loc(xr), T.get_indexer(pl.target[pl.regulator == xr])]
    assert np.mean(np.abs(tt) > 5) > 0.8


def test_homdel_indicator(sim, monkeypatch):
    """Deleted targets collapse (non-linear in x_B); a partner segment is amplified in the deleted lines.
    With the deletion indicator in W_B there are no partner -> deleted-target false positives; without it there are."""
    ds, tr, m1, b1 = sim
    R, T, _ = _index(ds)
    rows = R.isin(tr["homdel_partner_genes"])
    cols = T.isin(tr["homdel_genes"])
    assert rows.sum() >= 3 and cols.sum() >= 5
    assert _fp(m1["t"], rows, cols) < 0.05
    _ablate(monkeypatch, "homdel")
    assert _fp(run(ds, Params(), log=quiet)["t"], rows, cols) > 0.3


def test_cn_only_genes_and_a3_without_nan(sim):
    """CN-only genes (no expression column): t_cis = 0, s = +1, and A3 with w_co > 0 (+G2) adds no NaN
    (PI decision D5); G2 is a multiplier for every regulator, including singletons."""
    ds, tr, m1, _ = sim
    R, T, _ = _index(ds)
    cn_only = R.isin(tr["cn_only"])
    assert cn_only.sum() > 10
    assert (m1["cis_t"][cn_only] == 0).all() and (m1["s_cis"][cn_only] == 1).all()
    o = run(ds, Params(a3_rho=0.95, a3_wco=0.5, g2_threshold=3.0, g2_w=0.0), log=quiet)
    assert (np.isnan(o["score"]) == np.isnan(o["t"])).all()
    expressed = ds.context_expr["HCT116"].reindex(list(R)).to_numpy() >= 3.0
    assert (~expressed).sum() > 50
    assert (np.nan_to_num(o["score"][~expressed]) == 0).all()
    np.testing.assert_allclose(o["score"][expressed], (o["t"] * o["alpha"])[expressed], rtol=1e-5, equal_nan=True)


def test_hinge_x_uses_sex_centred_cn():
    """Hinge for X regulators splits at the sex-group median (copy-neutral), checked against statsmodels."""
    import statsmodels.api as sm
    from depmaptr.covariates import sex_centre, sex_groups
    from depmaptr.pipeline import GroupCache, group_t, prepare
    ds, tr = simulate(seed=6)
    prep = prepare(ds, Params(dosage="hinge"))
    xr = [k for k in ds.regulators if ds.genes.loc[k, "chrom"] == "X"][:3]
    gv = GroupCache(prep).get("X")
    t, df, _ = group_t(prep, gv, xr)
    D, _ = prep.builder.design("X")
    xs = sex_centre(prep.x_all[:, [prep.cn_idx[k] for k in xr]], sex_groups(ds.panel))
    for i in range(3):
        for j in [0, 100, 300]:
            cols = [D, np.maximum(xs[:, i], 0), prep.XT[:, j]] + ([prep.DEL[:, j]] if prep.DEL[:, j].any() else [])
            f = sm.OLS(prep.YT[:, j], np.column_stack(cols + [np.minimum(xs[:, i], 0)])).fit()
            assert abs(f.tvalues[-1] - t[i, j]) < 1e-8 * max(1, abs(f.tvalues[-1]))


def test_negative_cis_sign_rule(sim):
    ds, tr, m1, _ = sim
    R, T, _ = _index(ds)
    i = R.get_loc(tr["neg_cis_regulator"])
    assert m1["cis_t"][i] <= -3 and m1["s_cis"][i] == -1
    assert (m1["s_cis"][R.isin(tr["planted_regulators"]) & ~R.isin([tr["neg_cis_regulator"]])] == 1).all()
    pl = tr["planted"]
    ri, ti = R.get_indexer(pl.regulator), T.get_indexer(pl.target)
    pred = -np.sign(m1["score"][ri, ti]) * m1["s_cis"][ri]            # INTERFACES.md knockdown sign
    assert np.mean(pred == pl.kd_sign) > 0.95
    neg = (pl.regulator == tr["neg_cis_regulator"]).to_numpy()
    assert neg.sum() > 0 and np.mean(pred[neg] == pl.kd_sign[neg]) > 0.95
    # ignoring s_A would get the negative-cis regulator's signs wrong
    assert np.mean((-np.sign(m1["score"][ri, ti]))[neg] == pl.kd_sign[neg]) < 0.05


def test_cosegment_attribution(sim):
    ds, tr, _, _ = sim
    R, T, _ = _index(ds)
    a3 = run(ds, Params(a3_rho=0.95, a3_s=2.0), log=quiet)
    al = a3["alpha"]
    wins, chance = [], []
    pl = tr["planted"]
    for A, seg in tr["planted_segments"].items():
        members = [g for g in seg if g in R]
        ai = R.get_indexer(members)
        for b in pl.target[pl.regulator == A]:
            col = al[ai, T.get_loc(b)]
            wins.append(members[int(np.argmax(col))] == A)
            chance.append(1 / len(members))
    assert np.mean(wins) > 2 * np.mean(chance) and np.mean(wins) > 0.6
    # alpha sums to <= 1 within any segment and score = t * alpha
    np.testing.assert_allclose(a3["score"], a3["t"] * al, rtol=1e-5, equal_nan=True)


@pytest.mark.parametrize("seed", SEEDS)
def test_null_calibration(seed):
    ds, _ = simulate(seed=seed, null=True)
    R, T, trans = _index(ds)
    lam = lambda_gc(run(ds, Params(), log=quiet)["t"][trans])
    assert 0.9 < lam < 1.12, lam


@pytest.mark.parametrize("seed", SEEDS)
def test_nc1_calibration(seed):
    ds, _ = simulate(seed=seed)
    R, T, trans = _index(ds)
    lam = lambda_gc(run(ds, Params(nc1_index=seed), log=quiet)["t"][trans])
    assert 0.85 < lam < 1.12, lam


def test_variants_run_on_sim():
    """Every pre-registered model variant runs and returns finite trans scores of the right shape."""
    ds, tr = simulate(seed=4)
    R, T, trans = _index(ds)
    for kw in [dict(arm_handling="focal"), dict(arm_handling="loao"), dict(arm_handling="pcs", n_pcs=5),
               dict(g1_tau=10, g1_w=0.5), dict(g2_threshold=1.0, g2_w=0.0), dict(a3_rho=0.98, a3_s=5, a3_wco=0.5),
               dict(target_std=True), dict(dosage="hinge"), dict(cn_floor=2 ** -5), dict(rint=True)]:
        o = run(ds, Params(**kw), log=quiet)
        assert o["t"].shape == o["score"].shape == (len(R), len(T))
        assert np.isfinite(o["score"][trans]).mean() > 0.999, kw
        if kw.get("dosage") == "hinge":
            assert "t_gain" in o and np.isfinite(o["t_gain"][trans]).mean() > 0.99
        if kw.get("target_std"):
            z = o["score"][trans.copy()]
            assert abs(np.nanmean(z)) < 0.05
    for b in ["B0", "B1", "B2", "B3", "B4", "B4e", "B5", "B6"]:
        o = run_baseline(ds, b, log=quiet)
        assert o["score"].shape == (len(R), len(T)) and o["s_cis"].shape == (len(R),)


def test_focal_removes_arm_level_signal():
    """Focal adjusts A's own arm: the 3p driver's arm-wide LOCO signal mostly vanishes for non-driver 3p genes."""
    ds, tr = simulate(seed=5)
    R, T, _ = _index(ds)
    m = (ds.genes.loc[R, "arm_id"].to_numpy() == tr["driver_arm"]) & (R != tr["driver"])
    cols = T.isin(tr["driver_targets"])
    loco = run(ds, Params(), log=quiet)["t"][np.ix_(m, cols)]
    focal = run(ds, Params(arm_handling="focal"), log=quiet)["t"][np.ix_(m, cols)]
    assert np.mean(np.abs(loco) > 3) > 0.5 and np.mean(np.abs(focal) > 3) < np.mean(np.abs(loco) > 3) / 2
