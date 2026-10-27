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


@pytest.mark.parametrize("variant", ["M1", "B5", "a3_r095_s2_w0+g2_thr1_w0", "NC1"])
def test_run_model_sim(tmp_path, variant):
    out = tmp_path / variant.replace("+", "_")
    subprocess.run([sys.executable, str(ROOT / "scripts" / "run_model.py"), "--variant", variant, "--sim", "1",
                    "--out", str(out)], check=True, capture_output=True)
    z = _check_npz(out / "scores.npz")
    meta = json.loads((out / "meta.json").read_text())
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


@pytest.mark.parametrize("run_id", ["M1", "B0", "B1", "B2", "B3", "B4", "B4e", "B5", "B6", "NC1"])
def test_real_score_files(run_id):
    d = SCORES / run_id
    if not (d / "scores.npz").exists():
        pytest.skip(f"no real score file for {run_id}")
    meta = json.loads((d / "meta.json").read_text())
    assert META_KEYS <= set(meta)
    assert meta["n_samples"] == 1436
    assert 11000 < meta["n_regulators"] < 13000 and 17500 < meta["n_targets"] < 18500
    z = _check_npz(d / "scores.npz", meta["n_regulators"], meta["n_targets"])
    assert np.isfinite(z["score"]).mean() > 0.95   # B2: CN-only regulators are NaN
