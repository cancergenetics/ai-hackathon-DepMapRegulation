#!/usr/bin/env python
"""Collect runtimes from every score run's meta.json into results/runtime_table.tsv (PI decision D8).

    .venv/bin/python scripts/runtime_table.py [--scores data/results/scores] [--out results/runtime_table.tsv]

Columns: run_id, variant, kind, data, sample_set, n_samples, n_regulators, n_targets, wall_time_s (end of
scoring incl. the in-process NC1 call-set pass, excluding file writing), engine_s (pipeline.run of the main pass,
when recorded), calls_nc1_perms, cpu_count, blas_threads, peak_rss_gb (macOS ru_maxrss), git.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def row(meta: dict, run_dir: Path) -> dict:
    d = meta.get("diagnostics", {}) or {}
    th = meta.get("threads", {}) or {}
    blas = th.get("blas") or []
    calls = meta.get("calls") or {}
    return {
        "run_id": meta.get("run_id", run_dir.name), "variant": meta.get("variant"), "kind": meta.get("kind"),
        "data": meta.get("data"), "sample_set": meta.get("sample_set"), "n_samples": meta.get("n_samples"),
        "n_regulators": meta.get("n_regulators"), "n_targets": meta.get("n_targets"),
        "wall_time_s": meta.get("wall_time_s"),
        "engine_s": round(d["wall_engine_s"], 2) if isinstance(d.get("wall_engine_s"), (int, float)) else None,
        "calls_nc1_perms": calls.get("n_perm"), "cpu_count": th.get("cpu_count"),
        "blas_threads": ";".join(str(b.get("num_threads")) for b in blas) if blas else "n/a (Accelerate)",
        "peak_rss_gb": meta.get("peak_rss_gb"), "git": meta.get("git"),
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default=str(ROOT / "data" / "results" / "scores"))
    ap.add_argument("--out", default=str(ROOT / "results" / "runtime_table.tsv"))
    a = ap.parse_args(argv)
    rows = []
    for m in sorted(Path(a.scores).glob("*/meta.json")):
        try:
            rows.append(row(json.loads(m.read_text()), m.parent))
        except (json.JSONDecodeError, OSError) as e:
            print(f"skip {m}: {e}")
    df = pd.DataFrame(rows).sort_values(["kind", "run_id"]) if rows else pd.DataFrame()
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"{len(df)} runs -> {out}")


if __name__ == "__main__":
    main()
