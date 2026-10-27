#!/usr/bin/env python
"""The ONE final held-out (TEST) evaluation (DESIGN §6, §8.4). Run by the PI, once.

Guards (real run):
  1. results/final/selected.yaml exists, is tracked by git and has no uncommitted changes;
  2. results/final/LOCK exists neither in the working tree nor in HEAD (any rerun refuses);
  3. the sealed TEST files match the sha256 recorded in results/FROZEN.sha256.
Then it restores read access to the sealed TEST label files, evaluates every run listed in selected.yaml on
TEST (both lines), runs the paired comparisons selected-vs-each, removes read access again, writes
results/final/LOCK and commits results/final/.

--mock (the dry run DESIGN §6 requires before the real run): TRAIN fold-2 arms are relabelled as a
pseudo-TEST split and evaluated through the same code path, from the TRAIN labels only. Nothing is unsealed,
no LOCK is written, nothing is committed; output goes to results/final_mock/.

selected.yaml (all run ids refer to data/results/scores/<run_id>/):
  selected: M1_best            # the selected method
  ablations: [M1_best_noG2, M1_best_noG1, M1_best_focal]   # <= 3 (DESIGN §8.4)
  baselines: [B0, B1, B2, B3, B4, B4e, B5, B6]
  controls: [NC1]
  sensitivity: [M1_best_S1, ...]
  subsets:                     # optional regulator-subset evaluations, e.g. NC2
    - {name: NC2, run: M1_best_noG2, own_class: iii, expr_context_lt: 1.0}
  lines: [hct116, hek293t]     # optional; default both
  n_perm: 1000                 # optional
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import yaml

from depmaptr import benchmark as bm

ROOT = bm.ROOT
FINAL = ROOT / "results" / "final"
SELECTED = FINAL / "selected.yaml"
LOCK = FINAL / "LOCK"
FROZEN = ROOT / "results" / "FROZEN.sha256"
DEFAULT_SEALED = Path(os.environ.get("DEPMAPTR_SEALED", Path.home() / "DepMapTR_sealed"))


def git(*args, check=True) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=check)


def refuse(msg: str) -> None:
    raise SystemExit(f"final_eval REFUSED: {msg}")


def check_guards(mock: bool) -> None:
    if LOCK.exists() or git("cat-file", "-e", "HEAD:results/final/LOCK", check=False).returncode == 0:
        if not mock:
            refuse("results/final/LOCK is present: the TEST set has already been evaluated (DESIGN §6)")
    if not SELECTED.exists():
        refuse(f"{SELECTED.relative_to(ROOT)} does not exist")
    tracked = git("ls-files", "--error-unmatch", "results/final/selected.yaml", check=False).returncode == 0
    clean = git("diff", "--quiet", "HEAD", "--", "results/final/selected.yaml", check=False).returncode == 0
    if not (tracked and clean):
        if mock:
            print("WARNING (mock): selected.yaml is not committed; the real run will refuse", file=sys.stderr)
        else:
            refuse("results/final/selected.yaml must be committed (tracked, no local changes)")


def frozen_entries() -> dict:
    out = {}
    if FROZEN.exists():
        for ln in FROZEN.read_text().splitlines():
            if ln.strip() and not ln.startswith("#"):
                h, p = ln.split(None, 1)
                out[p.strip()] = h
    return out


def label_files(sealed: Path, lines) -> list[Path]:
    bd = sealed / "data" / "processed" / "bench"
    return [bd / f"{kind}_{line}_TEST.{ext}" for line in lines
            for kind, ext in (("positives", "parquet"), ("regulator_labels", "tsv"))]


def set_readable(paths, readable: bool) -> None:
    for p in paths:
        mode = p.stat().st_mode
        p.chmod(mode | stat.S_IRUSR if readable else mode & ~(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mock", action="store_true", help="dry run on TRAIN fold-2 arms relabelled as pseudo-TEST")
    ap.add_argument("--sealed-dir", type=Path, default=DEFAULT_SEALED)
    a = ap.parse_args()
    check_guards(a.mock)
    sel = yaml.safe_load(SELECTED.read_text())
    lines = sel.get("lines", list(bm.LINES))
    n_perm = int(sel.get("n_perm", bm.N_PERM))
    ablations = sel.get("ablations", []) or []
    if len(ablations) > 3:
        refuse("at most 3 ablations (DESIGN §8.4)")
    runs = [sel["selected"], *ablations, *(sel.get("baselines") or []), *(sel.get("controls") or []),
            *(sel.get("sensitivity") or [])]
    for r in runs + [s["run"] for s in sel.get("subsets") or []]:
        bm.resolve_run_dir(r)                           # fail before unsealing anything

    labels = label_files(a.sealed_dir, lines)
    if not a.mock:
        fe = frozen_entries()
        for p in labels:
            if not p.exists():
                refuse(f"sealed file missing: {p}")
            key = next((k for k in fe if Path(k).name == p.name), None)
            if key is None:
                refuse(f"{p.name} has no entry in results/FROZEN.sha256")
            set_readable([p], True)
            if bm.sha256_file(p) != fe[key]:
                set_readable([p], False)
                refuse(f"sha256 mismatch for {p}")
    token = bm.final_eval_token(a.sealed_dir / "data" / "processed" / "bench", mock=a.mock)
    out_dir = ROOT / "results" / ("final_mock" if a.mock else "final")
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = {"started": dt.datetime.now().isoformat(timespec="seconds"), "mock": a.mock,
               "git_hash": bm.git_hash(ROOT), "selected_yaml_sha256": bm.sha256_file(SELECTED), "results": [],
               "comparisons": []}
    try:
        for line in lines:
            for run in runs:
                out = bm.evaluate(run, line, "TEST", n_perm=n_perm, _token=token, _mock=a.mock)
                summary["results"].append({"run": run, "line": line, "split": out["split"],
                                           "primary_or": out["metrics"]["primary_or"],
                                           "perm_p": out["metrics"]["perm_regulator"]["perm_p"],
                                           "jackknife_ci95_log_or": out["metrics"]["jackknife_log_or"].get("ci95")})
            for s in sel.get("subsets") or []:
                regs = bm.load_bench(line, "TEST", token=token, mock=a.mock).regs
                mask = (regs["own_class"].astype(str) == str(s["own_class"])).values if "own_class" in s else True
                if "expr_context_lt" in s:
                    mask = mask & (regs["expr_context"].fillna(0).values < float(s["expr_context_lt"]))
                out = bm.evaluate(s["run"], line, "TEST", n_perm=n_perm, reg_subset=mask, subset_name=s["name"],
                                  _token=token, _mock=a.mock)
                summary["results"].append({"run": s["run"], "subset": s["name"], "line": line,
                                           "primary_or": out["metrics"]["primary_or"]})
            for other in runs[1:]:
                summary["comparisons"].append(bm.compare(sel["selected"], other, line, "TEST", _token=token,
                                                         _mock=a.mock))
    finally:
        if not a.mock:
            set_readable([p for p in labels if p.exists()], False)
    summary["finished"] = dt.datetime.now().isoformat(timespec="seconds")
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    if a.mock:
        print(f"mock final evaluation written to {out_dir}", file=sys.stderr)
        return 0
    LOCK.write_text(json.dumps({"locked": summary["finished"], "git_hash": summary["git_hash"],
                                "selected_yaml_sha256": summary["selected_yaml_sha256"]}, indent=1))
    git("add", "results/final")
    git("commit", "-m", "Final TEST evaluation (scripts/final_eval.py); results/final/LOCK written")
    print("final evaluation done; results/final/LOCK committed", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
