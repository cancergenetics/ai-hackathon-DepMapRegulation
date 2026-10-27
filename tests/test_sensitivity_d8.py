"""D8 model-side sensitivity analyses: hidden factors, CN-outlier exclusion, TP53 dropped for chr17."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from depmaptr import io as dio
from depmaptr.assoc import Residualiser
from depmaptr.pipeline import GroupCache, Params, hidden_factors, prepare, run
from depmaptr.simulate import simulate

ROOT = Path(__file__).resolve().parents[1]
SEL = dict(g1_tau=100, g1_w=0.0, a3_rho=0.95)
quiet = lambda *a, **k: None  # noqa: E731


def _run_model_module():
    spec = importlib.util.spec_from_file_location("run_model", ROOT / "scripts" / "run_model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_sensitivity_variants_compose_with_selected():
    rm = _run_model_module()
    defaults, table = rm.load_variants(ROOT / "results" / "variants.yaml")
    exp = {"S1": {"sample_set": "S1"}, "noHCT116": {"sample_set": "S2_noHCT116"},
           "hidden_factors": {"hidden_factors": 10}, "cn_outlier_exclusion": {"sample_set": "S2_noCNoutlier"},
           "tp53_chr17_drop": {"tp53_chr17_drop": True}}
    for vid, kv in exp.items():
        kind, params, _ = rm.resolve(f"g1_t100_w0+a3_r095_s2_w0+{vid}", defaults, table)
        assert kind == "model" and params["g1_tau"] == 100 and params["a3_rho"] == 0.95
        assert all(params[k] == v for k, v in kv.items())
        Params.from_dict(params)


# ----------------------------------------------------------------------------- hidden factors
def test_hidden_factors_orthogonal_to_design_and_chr_cn():
    ds, tr = simulate(seed=1)
    prep = prepare(ds, Params(hidden_factors=10))
    for c in ("2", "X"):
        D, _ = prep.builder.design(c)
        F, info = hidden_factors(prep, c, D)
        assert F.shape == (ds.n, 10) and info["hf_n_segments"] > 5
        keys = [k for k in ds.cn_keys if ds.genes.loc[k, "chrom"] == c]
        xc = prep.x_all[:, [prep.cn_idx[k] for k in keys]]
        Dl, _ = prep.builder.design(c)
        # F is exactly orthogonal to D_c; chr-c gene CN is (almost) fully explained by its segment means,
        # so F carries (almost) no chr-c CN: |corr| small for every chr-c gene
        assert np.abs(D.T @ F).max() < 1e-8 * np.abs(F).max() * ds.n
        Xr = Residualiser(Dl)(xc)
        corr = (F / np.linalg.norm(F, axis=0)).T @ (Xr / np.linalg.norm(Xr, axis=0))
        assert np.abs(corr).max() < 0.25
        # the group design carries the factors
        D2, names, extra = GroupCache(prep).design(c)
        assert D2.shape[1] == D.shape[1] + 10 and names[-1] == "hf10" and "hf_var_explained" in extra


def test_hidden_factors_keep_planted_signal(sim1=None):
    ds, tr = simulate(seed=2)
    R, T = pd.Index(ds.regulators), pd.Index(ds.targets)
    o0 = run(ds, Params(**SEL), log=quiet)
    o = run(ds, Params(**SEL, hidden_factors=10), log=quiet)
    pl = tr["planted"]
    ri, ti = R.get_indexer(pl.regulator), T.get_indexer(pl.target)
    assert np.mean(np.abs(o["t"][ri, ti]) > 5) > 0.9
    assert (np.isnan(o["score"]) == np.isnan(o["t"])).all()
    a, b = o0["t"].ravel(), o["t"].ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    assert np.corrcoef(a[ok], b[ok])[0, 1] > 0.9
    with pytest.raises(ValueError):
        run(ds, Params(hidden_factors=10, arm_handling="focal"), log=quiet)


# ----------------------------------------------------------------------------- TP53 chr17
def test_tp53_dropped_for_chr17_only():
    ds, tr = simulate(seed=3, n=300, n_auto=17, genes_per_arm=25, n_planted=6, targets_per_reg=10)
    p = Params(tp53_chr17_drop=True)
    prep = prepare(ds, p)
    gc = GroupCache(prep)
    _, n17, ex17 = gc.design("17")
    _, n16, _ = gc.design("16")
    assert "Z:tp53_altered" not in n17 and "Z:tp53_altered" in n16 and ex17["dropped"] == ["Z:tp53_altered"]
    o0 = run(ds, Params(), log=quiet)
    o1 = run(ds, p, log=quiet)
    rch = ds.genes.loc[ds.regulators, "chrom"].to_numpy()
    np.testing.assert_array_equal(o0["t"][rch != "17"], o1["t"][rch != "17"])
    assert not np.allclose(np.nan_to_num(o0["t"][rch == "17"]), np.nan_to_num(o1["t"][rch == "17"]))


# ----------------------------------------------------------------------------- CN-outlier exclusion
def test_cn_outlier_exclusion_sim():
    ds, tr = simulate(seed=4)
    sc = dio.cn_outlier_score(ds)
    # homozygous deletions alone do not count (floor 2^-2 -> |log2| <= 2)
    assert np.all(sc >= 0) and np.isfinite(sc).all()
    d2 = dio.apply_sample_set(ds, "S2_noCNoutlier")
    k = int(np.ceil(0.01 * ds.n))
    assert d2.n == ds.n - k and d2.sample_notes["cn_outlier_n_removed"] == k
    removed = {r["ModelID"] for r in d2.sample_notes["cn_outlier_removed"]}
    assert set(ds.model_ids[np.argsort(-sc)[:k]]) == removed and not removed & set(d2.model_ids)
    assert min(r["max_abs_log2_relcn"] for r in d2.sample_notes["cn_outlier_removed"]) >= np.sort(sc)[::-1][k]
    o = run(d2, Params(**SEL), log=quiet)
    assert o["t"].shape == (len(ds.regulators), len(ds.targets))


@pytest.mark.skipif(not (ROOT / "data" / "processed" / "panel_s2.tsv").exists(), reason="real panel absent")
def test_cn_outlier_exclusion_real():
    ds = dio.load_real("S2_noCNoutlier", checksums=False)
    notes = ds.sample_notes
    assert notes["cn_outlier_n_removed"] == 15 and ds.n == 1421
    names = {r["name"] for r in notes["cn_outlier_removed"]}
    assert "COLO 320" in names or any("COLO" in n and "320" in n for n in names), names
    assert notes["cn_outlier_threshold"] > 2          # amplicons, not deletions (floor gives at most 2)
