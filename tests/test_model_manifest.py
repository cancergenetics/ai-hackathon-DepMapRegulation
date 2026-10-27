"""Score-file contract (docs/INTERFACES.md) and manifest shapes (DESIGN §8.7)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SCORES = ROOT / "data" / "results" / "scores"
META_KEYS = {"run_id", "variant", "params", "git", "input_checksums", "n_samples", "wall_time_s", "threads",
             "peak_rss_gb"}


def _check_npz(path: Path, R=None, T=None):
    z = np.load(path)
    for k in ("regulators", "targets", "t", "score", "s_cis", "cis_t"):
        assert k in z.files, k
    R = R or len(z["regulators"])
    T = T or len(z["targets"])
    assert z["regulators"].dtype.kind == "U" and z["targets"].dtype.kind == "U"
    assert len(set(z["regulators"])) == R and len(set(z["targets"])) == T
    assert z["t"].shape == z["score"].shape == (R, T)
    assert z["t"].dtype == np.float32 and z["score"].dtype == np.float32
    assert z["s_cis"].shape == (R,) and z["s_cis"].dtype == np.int8 and set(np.unique(z["s_cis"])) <= {-1, 1}
    assert z["cis_t"].shape == (R,) and z["cis_t"].dtype == np.float32
    return z


@pytest.mark.parametrize("variant", ["M1", "B0", "B5", "a3_r095_s2_w0+g2_thr1_w0", "NC1"])
def test_run_model_sim(tmp_path, variant):
    out = tmp_path / variant.replace("+", "_")
    subprocess.run([sys.executable, str(ROOT / "scripts" / "run_model.py"), "--variant", variant, "--sim", "1",
                    "--out", str(out)], check=True, capture_output=True)
    z = _check_npz(out / "scores.npz")
    meta = json.loads((out / "meta.json").read_text())
    if variant in ("M1", "a3_r095_s2_w0+g2_thr1_w0"):
        assert "call" in z.files and z["call"].dtype == bool and z["call"].shape == z["t"].shape
        assert meta["calls"]["n_calls"] == int(z["call"].sum())
        rch = np.array([r.split("(")[0][1:-5] for r in z["regulators"]])
        tch = np.array([r.split("(")[0][1:-5] for r in z["targets"]])
        assert not z["call"][rch[:, None] == tch[None, :]].any()           # trans only
        assert np.all(np.abs(z["t"][z["call"]]) >= meta["calls"]["t_star"])
    else:
        assert "call" not in z.files
    if variant == "B0":
        assert meta["seed"] == 20261027
    assert META_KEYS <= set(meta) and meta["variant"] == variant and meta["data"] == "SIM"
    man = json.loads((out / "manifest.json").read_text())
    assert man["n_regulators"] == len(z["regulators"]) and man["n_targets"] == len(z["targets"])
    assert man["n_pairs"] == z["t"].size


def test_variants_yaml():
    v = yaml.safe_load((ROOT / "results" / "variants.yaml").read_text())
    ids = [e["id"] for e in v["variants"]]
    assert len(ids) == len(set(ids)) and len(ids) <= v["cap_train_configurations"]
    tiers = {e["tier"] for e in v["variants"]}
    assert tiers == set(range(8))
    assert {e["id"] for e in v["baselines"]} == {"B0", "B1", "B2", "B3", "B4", "B4e", "B5", "B6"}
    from depmaptr.pipeline import Params
    Params.from_dict(v["defaults"])
    for e in v["variants"]:
        Params.from_dict({**v["defaults"], **e["params"]})


RUN_DIRS = sorted(p.name for p in SCORES.glob("*") if (p / "scores.npz").exists()) if SCORES.exists() else []


@pytest.mark.parametrize("run_id", RUN_DIRS or ["<none>"])
def test_real_score_files(run_id):
    """Every real run directory: contract, pinned ranges, and NaN share <= A=B (+ CN-only rows for B2/B6)."""
    d = SCORES / run_id
    if not (d / "scores.npz").exists():
        pytest.skip("no real score files")
    meta = json.loads((d / "meta.json").read_text())
    assert META_KEYS <= set(meta)
    if meta.get("data") != "real":
        pytest.skip("not a real-data run")
    assert meta["n_samples"] == 1436
    assert 11000 < meta["n_regulators"] < 13000 and 17500 < meta["n_targets"] < 18500
    z = _check_npz(d / "scores.npz", meta["n_regulators"], meta["n_targets"])
    R, T = z["t"].shape
    n_same = len(set(z["regulators"]) & set(z["targets"]))
    allowed = n_same
    if meta["variant"] == "B2":       # whole NaN rows: CN-only regulators + constant-expression regulators
        nan_rows = np.isnan(z["score"]).all(1)
        assert int(np.isnan(z["score"]).sum()) == int(nan_rows.sum()) * T
        assert nan_rows.sum() <= meta["n_regulators_without_expression"] + 5
        allowed = int(nan_rows.sum()) * T
    if meta["variant"] == "B6":
        allowed = int(np.isnan(z["score"]).any(1).sum()) * T      # regulators on arms without an arm covariate
    assert int(np.isnan(z["score"]).sum()) <= allowed + 100, (int(np.isnan(z["score"]).sum()), allowed)
    if meta["kind"] == "model" and meta["nc1_index"] is None:
        assert "call" in z.files and meta["calls"]["n_calls"] == int(z["call"].sum())
