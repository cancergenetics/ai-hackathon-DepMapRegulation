"""Model-side invariants (DESIGN §8.6) and unit tests for covariates, gates and attribution."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from depmaptr import gates
from depmaptr.attribution import greedy_segments, log_abf, segment_alpha
from depmaptr.covariates import build_Z, lineage_levels

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
real = pytest.mark.skipif(not (PROC / "panel_s2.tsv").exists(), reason="real panel absent")


# ----------------------------------------------------------------------------- gates / sign rule
def test_s_cis_rule_and_kd_sign():
    t = np.array([-10.0, -3.0, -2.99, 0.0, 5.0, np.nan])
    s = gates.s_cis(t)
    assert s.tolist() == [-1, -1, 1, 1, 1, 1]
    # negative-cis regulator: beta_AB < 0 and s_A = -1 -> knockdown predicted to DECREASE B
    score = np.array([[-2.0], [-2.0], [2.0]])
    kd = gates.predicted_kd_sign(score, np.array([-1, 1, 1]))
    assert kd.ravel().tolist() == [-1.0, 1.0, -1.0]


def test_g1_g2():
    assert gates.g1_pass(np.array([4.0, 3.0, -5.0, 0.0]), 10).tolist() == [True, False, False, False]
    assert gates.g1_weight(np.array([4.0, 1.0]), 10, 0.5).tolist() == [1.0, 0.5]
    ctx = pd.Series({"A": 3.0, "B": 0.5})
    assert gates.g2_weight(["A", "B", "C"], ctx, 1.0, 0.0).tolist() == [1.0, 0.0, 0.0]   # C: no expression column


# ----------------------------------------------------------------------------- attribution
def test_greedy_segments_known_answer():
    rng = np.random.default_rng(0)
    n = 400
    a, b = rng.normal(size=n), rng.normal(size=n)
    X = np.column_stack([a + 0.01 * rng.normal(size=n) for _ in range(4)] +
                        [b + 0.01 * rng.normal(size=n) for _ in range(3)] + [rng.normal(size=n)])
    assert greedy_segments(X, 0.95).tolist() == [0, 0, 0, 0, 1, 1, 1, 2]


def test_segment_alpha_formula():
    labels = np.array([0, 0, 0, 1])
    t = np.array([[5.0, 0.0], [3.0, 0.0], [1.0, 0.0], [2.0, 2.0]])
    df = np.full_like(t, 1000.0)
    tcis = np.array([10.0, 0.0, 0.0, 1.0])
    al = segment_alpha(labels, t, df, tcis, s=2.0)
    lp = tcis[:3] / 2 + log_abf(t[:3, 0], 1000.0)
    ref = np.exp(lp - lp.max()) / np.exp(lp - lp.max()).sum()
    np.testing.assert_allclose(al[:3, 0], ref, rtol=1e-12)
    np.testing.assert_allclose(al[:3].sum(0), 1.0)
    assert (al[3] == 1).all()                                   # singleton
    # all pi = 0 (w_G2 = 0 for every member) -> uniform prior
    al0 = segment_alpha(labels, t, df, tcis, s=2.0, w_g2=np.array([0.0, 0.0, 0.0, 1.0]))
    lp0 = log_abf(t[:3, 0], 1000.0)
    np.testing.assert_allclose(al0[:3, 0], np.exp(lp0) / np.exp(lp0).sum(), rtol=1e-12)
    # pi_g with w_G2: a zero-weight member gets alpha 0
    al1 = segment_alpha(labels, t, df, tcis, s=2.0, w_g2=np.array([0.0, 1.0, 1.0, 1.0]))
    assert al1[0, 0] == 0 and abs(al1[1:3, 0].sum() - 1) < 1e-12


# ----------------------------------------------------------------------------- covariates
def _toy_panel():
    rows = []
    for L, d, k in [("Myeloid", "AML", 20), ("Myeloid", "MPN", 5), ("Lung", "NSCLC", 30), ("Thyroid", "T", 10),
                    ("Esophagus/Stomach", "EAC", 16), ("Esophagus/Stomach", "ESCC", 15), (None, None, 2)]:
        rows += [(L, d)] * k
    p = pd.DataFrame(rows, columns=["lineage", "primary_disease"])
    n = len(p)
    rng = np.random.default_rng(0)
    p["sex"] = rng.choice(["Male", "Female", "Unknown", None], n)
    p["ploidy"] = np.where(rng.random(n) < 0.3, np.nan, rng.normal(3, 0.5, n))
    p["ploidy_missing"] = p["ploidy"].isna()
    p["wgd"] = np.where(p["ploidy_missing"], np.nan, rng.integers(0, 2, n))
    p["wgd_missing"] = p["ploidy_missing"]
    p["msi_score"] = np.where(rng.random(n) < 0.2, np.nan, rng.gamma(1, 10, n))
    p["tp53_altered"] = np.where(p["msi_score"].isna(), None, rng.random(n) < 0.5)
    p["media_base"] = rng.choice(["RPMI", "DMEM", "X1", None], n, p=[0.5, 0.3, 0.1, 0.1])
    p["serum_free"] = np.where(p["media_base"].isna(), None, rng.random(n) < 0.1)
    p["growth_pattern"] = rng.choice(["Adherent", "Suspension"], n)
    p["stranded"] = rng.random(n) < 0.5
    p["platform"] = None
    return p


def test_lineage_rule():
    lv = lineage_levels(_toy_panel())
    vc = lv.value_counts()
    assert vc["Myeloid::AML"] == 20 and vc["Myeloid"] == 5        # remainder stays in the parent lineage
    assert vc["Esophagus/Stomach::EAC"] == 16 and vc["Esophagus/Stomach::ESCC"] == 15
    assert vc["other"] == 12                                     # Thyroid (10) + missing (2)
    assert vc["Lung"] == 30


def test_build_Z():
    Z, log = build_Z(_toy_panel())
    assert np.isfinite(Z.to_numpy()).all()
    assert "wgd_missing" in log["dropped_duplicate"]             # identical to ploidy_missing
    assert not any(c.startswith("platform") for c in Z.columns)  # PI decision D3
    assert {"msi_high", "msi_missing", "tp53_altered", "sex[Female]", "sex[Unknown]", "serum_free", "stranded"} <= set(Z.columns)
    assert not Z.T.duplicated().any()


# ----------------------------------------------------------------------------- real data
@pytest.fixture(scope="module")
def real_ds():
    from depmaptr import io as dio
    return dio.load_real("S2", checksums=False)


@real
def test_real_invariants(real_ds):
    from depmaptr import io as dio
    inv = dio.check_invariants(real_ds, real=True)
    assert inv["median_cis_r"] > 0.1
    assert inv["n_samples"] == 1436
    assert 11000 < inv["n_regulators"] < 13000 and 17500 < inv["n_targets"] < 18500
    assert len(real_ds.panel) == real_ds.relcn.shape[0] == real_ds.expr.shape[0]
    assert (real_ds.panel["ModelID"].to_numpy() == pd.read_csv(PROC / "panel_s2.tsv", sep="\t")["ModelID"].to_numpy()).all()
    # matrices are aligned to the panel: spot-check against the parquet
    cn = pd.read_parquet(PROC / "depmap_cn_matched.parquet", columns=["ModelID", real_ds.cn_keys[100]]).set_index("ModelID")
    np.testing.assert_allclose(real_ds.relcn[:, 100], cn.loc[real_ds.model_ids].iloc[:, 0].to_numpy(), rtol=1e-6)


@real
def test_real_context_lines(real_ds):
    from depmaptr import io as dio
    dio.context_model_check()
    assert "HCT116" in real_ds.context_expr
    ex = pd.read_parquet(PROC / "depmap_expr_matched.parquet", columns=["ModelID", "MYC (4609)"]).set_index("ModelID")
    assert real_ds.context_expr["HCT116"]["MYC (4609)"] == pytest.approx(ex.loc["ACH-000971", "MYC (4609)"])


@real
def test_real_arm_cn(real_ds):
    from depmaptr.covariates import arm_cn
    arms = arm_cn(real_ds)
    assert arms.shape == (1436, 41) and np.isfinite(arms.to_numpy()).all()
    # Autosomal arms agree with the data agent's derived parquet (same median-of-genes definition)
    der = pd.read_parquet(PROC / "depmap_arm_cn_derived.parquet")
    der = der.pivot(index="ModelID", columns="arm_id", values="median_log2_cn").loc[real_ds.model_ids]
    auto = [a for a in arms.columns if not a.startswith("X")]
    d = np.abs(arms[auto].to_numpy() - der[auto].to_numpy())
    assert np.median(d) < 1e-3 and np.percentile(d, 99) < 0.05
    # sex-centred X arms: median 0 within each sex group
    from depmaptr.covariates import sex_groups
    g = sex_groups(real_ds.panel)
    for s in np.unique(g):
        assert np.abs(np.median(arms.loc[g == s, ["Xp", "Xq"]].to_numpy(), axis=0)).max() < 1e-12


@real
def test_real_Z(real_ds):
    from depmaptr.assoc import Residualiser
    from depmaptr.covariates import DesignBuilder, arm_cn
    Z, log = build_Z(real_ds.panel)
    assert np.isfinite(Z.to_numpy()).all() and Z.shape[0] == 1436
    lv = log["lineage_levels"]
    assert lv["Myeloid::Acute Myeloid Leukemia"] == 40 and "Thyroid" not in lv
    b = DesignBuilder(real_ds, Z.to_numpy(), arm_cn(real_ds), "loco")
    for c in ["1", "17", "X"]:
        res = Residualiser(b.design(c)[0])          # asserts cond < 1e6
        assert res.rank == res.n_cols
