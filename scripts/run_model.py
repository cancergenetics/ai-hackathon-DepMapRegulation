#!/usr/bin/env python
"""Score all frozen variable-CN regulators x all measured targets for one variant (docs/INTERFACES.md).

    .venv/bin/python scripts/run_model.py --variant M1 --out data/results/scores/M1
    .venv/bin/python scripts/run_model.py --variant focal+g1_t30_w0 --out data/results/scores/focal_g1
    .venv/bin/python scripts/run_model.py --variant B1 --out data/results/scores/B1
    .venv/bin/python scripts/run_model.py --variant NC1 --perm-index 3 --out data/results/scores/NC1_3
    .venv/bin/python scripts/run_model.py --variant M1 --sim 1 --out /tmp/sim_M1     # SIM data

Writes <out>/scores.npz (regulators, targets, t, score, s_cis, cis_t [+ alpha, t_gain]), meta.json and
manifest.json (DESIGN §8.7).
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from depmaptr import io as dio  # noqa: E402
from depmaptr.baselines import run_baseline  # noqa: E402
from depmaptr import calls as dcalls  # noqa: E402
from depmaptr.pipeline import ROOT_SEED, Params, run  # noqa: E402

MEM_LIMIT_GB = 20.0


def load_variants(path: Path) -> tuple[dict, dict]:
    v = yaml.safe_load(path.read_text())
    table = {}
    for sec in ("variants", "baselines", "controls", "sensitivity", "posthoc_D6", "controls_D6"):
        for e in v.get(sec, []) or []:
            table[e["id"]] = e
    return v["defaults"], table


def resolve(spec: str, defaults: dict, table: dict) -> tuple[str, dict, str | None]:
    """'a+b' -> (kind, params, baseline id or None)."""
    params = dict(defaults)
    baseline = None
    for vid in spec.split("+"):
        if vid not in table:
            raise SystemExit(f"unknown variant id {vid!r} (see results/variants.yaml)")
        e = table[vid]
        if e.get("implemented") is False:
            raise SystemExit(f"variant {vid} is not implemented")
        if e.get("kind") == "baseline":
            baseline = vid
        params.update(e.get("params", {}) or {})
    return ("baseline" if baseline else "model"), params, baseline


def git_hash() -> str:
    try:
        h = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "src", "scripts"], cwd=ROOT, capture_output=True,
                               text=True).stdout.strip()
        return h + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def blas_info() -> dict:
    try:
        from threadpoolctl import threadpool_info
        return [{k: d.get(k) for k in ("internal_api", "num_threads", "version")} for d in threadpool_info()]
    except Exception:
        return []


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--perm-index", type=int, default=None, help="NC1 spawn index (overrides params)")
    ap.add_argument("--seed", type=int, default=None, help="B0 seed")
    ap.add_argument("--sim", type=int, default=None, help="run on SIM data with this seed instead of real data")
    ap.add_argument("--variants-file", default=str(ROOT / "results" / "variants.yaml"))
    ap.add_argument("--max-wall", type=int, default=3600, help="wall-clock limit in seconds (DESIGN §8.10)")
    ap.add_argument("--chunk", type=int, default=512)
    ap.add_argument("--nc1-perms", type=int, default=1,
                    help="model runs: NC1 passes (spawn indices 0..k-1) of the same variant for the call set")
    ap.add_argument("--no-calls", action="store_true", help="skip the call set (no in-process NC1)")
    a = ap.parse_args(argv)

    def _timeout(*_):
        raise TimeoutError(f"wall-clock limit {a.max_wall}s exceeded")
    signal.signal(signal.SIGALRM, _timeout)
    signal.alarm(a.max_wall)

    t0 = time.time()
    defaults, table = load_variants(Path(a.variants_file))
    kind, params, baseline = resolve(a.variant, defaults, table)
    if a.perm_index is not None:
        params["nc1_index"] = a.perm_index
    params["chunk"] = a.chunk
    p = Params.from_dict(params)

    if a.sim is not None:
        from depmaptr.simulate import simulate
        ds, _ = simulate(seed=a.sim)
        inv = dio.check_invariants(ds, real=False)
    else:
        ds = dio.load_real(p.sample_set)
        inv = dio.check_invariants(ds, real=True)
    R, T, n = len(ds.regulators), len(ds.targets), ds.n
    # compute guard: largest per-chunk working set (k+2)^2/2 R_chunk x T arrays + residualised targets
    est_gb = (p.chunk * T * 8 * 12 + n * T * 8 * 4 + R * T * 4 * 4) / 1e9
    if p.a3_rho is not None:
        est_gb += 2200 * T * 8 * 6 / 1e9
    print(f"[run_model] {a.variant}: n={n} R={R} T={T}; est. peak {est_gb:.1f} GB", flush=True)
    assert est_gb < MEM_LIMIT_GB, f"estimated memory {est_gb:.1f} GB >= {MEM_LIMIT_GB} GB"

    call_meta = None
    if kind == "baseline":
        seed_used = (ROOT_SEED if a.seed is None else a.seed) if baseline == "B0" else None
        out = run_baseline(ds, baseline, p, seed=seed_used)
    else:
        seed_used = None
        out = run(ds, p)
        if p.nc1_index is None and not a.no_calls:
            # MVP call set (DESIGN §12) against NC1 passes of the SAME variant. Gates, A3 and target
            # standardisation do not change t, so they are switched off in the null passes.
            gc = ds.genes
            trans = dcalls.trans_mask(gc.loc[ds.regulators, "chrom"].astype(str).to_numpy(),
                                      gc.loc[ds.targets, "chrom"].astype(str).to_numpy())
            real_h = dcalls.abs_t_hist(out["t"], trans)
            null_h = np.zeros_like(real_h)
            pn = Params.from_dict({**p.to_dict(), "g1_tau": None, "g2_threshold": None, "a3_rho": None,
                                   "target_std": False, "score_base": "t", "fuse_gamma": None, "regctx": False})
            for k in range(a.nc1_perms):
                pn.nc1_index = k
                null_h += dcalls.abs_t_hist(run(ds, pn, log=lambda *x, **y: None)["t"], trans)
            ts = dcalls.t_star(real_h, null_h, a.nc1_perms)
            out["call"] = dcalls.call_matrix(out["t"], trans, ts["t_star"])
            ts["n_calls"] = int(out["call"].sum())
            call_meta = {**ts, "scope": "trans pairs", "nc1_source": "in-process NC1 passes of this variant",
                         "nc1_indices": list(range(a.nc1_perms)), "nc1_params": pn.to_dict() | {"nc1_index": "0..k-1"},
                         "fdr": dcalls.FDR, "bin_width": dcalls.BIN}
            del trans
            print(f"[run_model] call set: t* = {ts['t_star']}, {ts['n_calls']} trans pairs", flush=True)
    wall = time.time() - t0

    od = Path(a.out)
    od.mkdir(parents=True, exist_ok=True)
    arrays = {"regulators": np.asarray(ds.regulators, dtype=str), "targets": np.asarray(ds.targets, dtype=str),
              "t": np.asarray(out["t"], np.float32), "score": np.asarray(out["score"], np.float32),
              "s_cis": np.asarray(out["s_cis"], np.int8), "cis_t": np.asarray(out["cis_t"], np.float32)}
    for extra in ("alpha", "t_gain", "pair_score", "target_weight"):
        if extra in out:
            arrays[extra] = np.asarray(out[extra], np.float32)
    if "call" in out:
        arrays["call"] = np.asarray(out["call"], bool)
    assert arrays["t"].shape == arrays["score"].shape == (R, T)
    assert arrays["s_cis"].shape == arrays["cis_t"].shape == (R,)
    np.savez(od / "scores.npz", **arrays)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_rss_gb = rss / 1e9 if sys.platform == "darwin" else rss * 1024 / 1e9
    diag = out.get("diagnostics", {})
    n_trans = int(sum(((ds.genes.loc[ds.regulators, "chrom"].astype(str).to_numpy() == c).sum()
                       * (ds.genes.loc[ds.targets, "chrom"].astype(str).to_numpy() != c).sum())
                      for c in dio.CHROMS))
    n_nan = int(np.isnan(arrays["score"]).sum())
    meta = {
        "run_id": od.name, "variant": a.variant, "kind": kind, "params": p.to_dict(),
        "seed": seed_used, "nc1_index": p.nc1_index, "calls": call_meta,
        "git": git_hash(), "data": "SIM" if a.sim is not None else "real", "sim_seed": a.sim,
        "input_checksums": ds.input_checksums, "sample_set": ds.sample_set,
        "n_samples": int(n), "n_regulators": int(R), "n_targets": int(T), "n_pairs": int(R * T),
        "n_trans_pairs": n_trans, "n_nan_scores": n_nan,
        "n_regulators_without_expression": int(sum(k not in set(ds.expr_keys) for k in ds.regulators)),
        "invariants": inv, "wall_time_s": round(wall, 2), "threads": {"cpu_count": os.cpu_count(), "blas": blas_info()},
        "peak_rss_gb": round(peak_rss_gb, 3), "python": platform.python_version(), "numpy": np.__version__,
        "diagnostics": diag,
    }
    (od / "meta.json").write_text(json.dumps(meta, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    manifest = {k: meta[k] for k in ("run_id", "variant", "n_samples", "n_regulators", "n_targets", "n_pairs",
                                     "n_trans_pairs", "sample_set")}
    manifest["n_models"] = manifest.pop("n_samples")
    manifest["covariates"] = diag.get("Z_columns", [])
    manifest["n_arms"] = diag.get("n_arms")
    (od / "manifest.json").write_text(json.dumps(manifest, indent=1))
    signal.alarm(0)
    print(f"[run_model] done {a.variant} in {wall:.1f}s, peak RSS {peak_rss_gb:.2f} GB -> {od}", flush=True)


if __name__ == "__main__":
    main()
