#!/usr/bin/env python
"""The ONE final held-out (TEST) evaluation (DESIGN §6, §8.4). Run by the PI, once.

Order of a real run:
  1. guards: results/final/LOCK absent (working tree and HEAD); results/final/selected.yaml tracked and
     unmodified; results/FROZEN.sha256 present with no drift in any in-tree frozen file; no TEST result
     JSON under results/final/ and no TEST row in results/tuning_log.tsv; <= 3 ablations; all runs exist;
  2. results/final/LOCK is written (status STARTED, timestamp, selected.yaml sha256, git hash) and
     COMMITTED, before any TEST file is touched: a crash after this point still counts as the one run;
  3. read access is restored on the sealed TEST label files only, and their sha256 is checked against the
     sealed block of results/FROZEN.sha256;
  4. every run in selected.yaml is evaluated on TEST (both lines), plus the regulator subsets (e.g. NC2),
     the split-sibling sensitivity of the selected run, PC2, and the paired comparisons selected-vs-each;
  5. read access is removed again; LOCK is updated to COMPLETED (or FAILED with the error) and
     results/final/ is committed.

--mock (the dry run DESIGN §6 requires): builds a mock sealed directory (data/interim/mock_sealed/) holding
TRAIN fold-2 labels under the real TEST file names, chmod a-rw, with its own FROZEN.sha256, and then runs
steps 1-5 through the same code path (token, sealed-file reads, hash checks, chmod round trip, LOCK), with
output in results/final_mock/ and no git commit. Real TEST labels are never touched.

selected.yaml (run ids refer to data/results/scores/<run_id>/):
  selected: M1_best
  ablations: [M1_best_noG2, M1_best_noG1, M1_best_focal]     # <= 3 (DESIGN §8.4)
  baselines: [B0, B1, B2, B3, B4, B4e, B5, B6]
  controls: [NC1]
  sensitivity: [M1_best_S1, ...]
  subsets:                                                   # optional; subset in depmaptr.benchmark.SUBSETS
    - {name: NC2, run: M1_best_noG2, subset: nc2}
  lines: [hct116, hek293t]                                   # optional; default both
  n_perm: 1000                                               # optional
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

from depmaptr import benchmark as bm

ROOT = bm.ROOT
SELECTED = ROOT / "results" / "final" / "selected.yaml"
FROZEN = ROOT / "results" / "FROZEN.sha256"
DEFAULT_SEALED = Path(os.environ.get("DEPMAPTR_SEALED", Path.home() / "DepMapTR_sealed"))
MOCK_SEALED = ROOT / "data" / "interim" / "mock_sealed"


def git(*args, check=True) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=check)


def refuse(msg: str) -> None:
    raise SystemExit(f"final_eval REFUSED: {msg}")


def now() -> str:
    return dt.datetime.now().isoformat(timespec="seconds")


def set_readable(paths, readable: bool) -> None:
    for p in paths:
        mode = p.stat().st_mode
        p.chmod(mode | stat.S_IRUSR if readable else mode & ~(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH))


class Ctx:
    """Real and mock runs differ only in these locations and in whether git commits happen."""

    def __init__(self, mock: bool, sealed_dir: Path):
        self.mock = mock
        self.out = ROOT / "results" / ("final_mock" if mock else "final")
        self.lock = self.out / "LOCK"
        self.sealed = MOCK_SEALED if mock else sealed_dir
        self.frozen = (MOCK_SEALED / "FROZEN.sha256") if mock else FROZEN

    def commit(self, msg: str) -> None:
        if self.mock:
            print(f"(mock) would commit results/final: {msg}", file=sys.stderr)
            return
        git("add", "results/final")
        git("commit", "-m", msg)


def check_guards(ctx: Ctx, sel: dict, lines: list) -> None:
    if not ctx.mock:
        if ctx.lock.exists() or git("cat-file", "-e", "HEAD:results/final/LOCK", check=False).returncode == 0:
            refuse("results/final/LOCK is present: the TEST set has already been evaluated (DESIGN §6)")
        tracked = git("ls-files", "--error-unmatch", "results/final/selected.yaml", check=False).returncode == 0
        clean = git("diff", "--quiet", "HEAD", "--", "results/final/selected.yaml", check=False).returncode == 0
        if not (tracked and clean):
            refuse("results/final/selected.yaml must be committed (tracked, no local changes)")
        if any(ctx.out.glob("*/*.json")) or (ctx.out / "summary.json").exists():
            refuse("TEST result files already exist under results/final/")
        log = ROOT / "results" / "tuning_log.tsv"
        if log.exists() and (pd.read_csv(log, sep="\t", usecols=["split"])["split"] == "TEST").any():
            refuse("results/tuning_log.tsv already holds TEST rows")
        if not FROZEN.exists():
            refuse("results/FROZEN.sha256 missing: seal first (scripts/benchmark/seal_test.sh)")
    if FROZEN.exists():
        drift = bm.frozen_drift(FROZEN, ROOT)
        if drift:
            if ctx.mock:
                print(f"WARNING (mock): frozen files drifted {drift}; the real run will refuse", file=sys.stderr)
            else:
                refuse(f"frozen in-tree files changed since the freeze (DESIGN §8.5): {drift}")
    elif ctx.mock:
        print("WARNING (mock): results/FROZEN.sha256 not found (not sealed yet)", file=sys.stderr)
    ab = sel.get("ablations") or []
    if len(ab) > 3:
        refuse("at most 3 ablations (DESIGN §8.4)")
    for s in sel.get("subsets") or []:
        if s.get("subset", str(s["name"]).lower()) not in bm.SUBSETS:
            refuse(f"unknown subset in selected.yaml: {s}")
    for r in all_runs(sel) + [s["run"] for s in sel.get("subsets") or []]:
        bm.resolve_run_dir(r)                                   # fail before the LOCK is written


def all_runs(sel: dict) -> list:
    return [sel["selected"], *(sel.get("ablations") or []), *(sel.get("baselines") or []),
            *(sel.get("controls") or []), *(sel.get("sensitivity") or [])]


def label_paths(ctx: Ctx, lines) -> list:
    return [ctx.sealed / "data" / "processed" / "bench" / f"{kind}_{line}_TEST.{ext}" for line in lines
            for kind, ext in (("positives", "parquet"), ("regulator_labels", "tsv"))]


def unseal_and_verify(ctx: Ctx, labels: list) -> None:
    _, sealed = bm.read_frozen(ctx.frozen)
    want = {q: h for h, q in sealed}
    for p in labels:
        rel = str(p.relative_to(ctx.sealed))
        if not p.exists():
            raise RuntimeError(f"sealed file missing: {p}")
        if rel not in want:
            raise RuntimeError(f"{rel} has no sealed entry in {ctx.frozen}")
        set_readable([p], True)
        if bm.sha256_file(p) != want[rel]:
            raise RuntimeError(f"sha256 mismatch for {p}")


def run_all(ctx: Ctx, sel: dict, lines: list, token, summary: dict) -> None:
    n_perm = int(sel.get("n_perm", bm.N_PERM))
    runs = all_runs(sel)
    for line in lines:
        for run in runs:
            out = bm.evaluate(run, line, "TEST", n_perm=n_perm, _token=token)
            summary["results"].append({"run": run, "line": line, "split": out["split"],
                                       "primary_or": out["metrics"]["primary_or"],
                                       "perm_p": out["metrics"]["perm_regulator"]["perm_p"],
                                       "jackknife_ci95_log_or": out["metrics"]["jackknife_log_or"].get("ci95")})
        regs = bm.load_bench(line, "TEST", token=token).regs
        subs = [dict(s, subset=s.get("subset", str(s["name"]).lower())) for s in sel.get("subsets") or []]
        subs.append({"name": "no_split_sibling", "run": sel["selected"], "subset": "no_split_sibling"})
        for s in subs:
            mask = bm.subset_mask(regs, s["subset"], ROOT)
            out = bm.evaluate(s["run"], line, "TEST", n_perm=n_perm, reg_subset=mask, subset_name=s["name"],
                              _token=token)
            summary["results"].append({"run": s["run"], "subset": s["name"], "line": line,
                                       "n_regulators": int(mask.sum()), "primary_or": out["metrics"]["primary_or"],
                                       "perm_p": out["metrics"]["perm_regulator"]["perm_p"]})
        for other in runs[1:]:
            summary["comparisons"].append(bm.compare(sel["selected"], other, line, "TEST", _token=token))
    if {"hct116", "hek293t"} <= set(lines):
        out = bm.evaluate_pc2("TEST", ROOT, n_perm=n_perm, _token=token)
        summary["results"].append({"run": "PC2", "line": "hct116", "primary_or": out["metrics"]["primary_or"],
                                   "callset_or": out["metrics"]["callset"].get("or")})
    else:
        summary["notes"].append("PC2 skipped: needs both lines")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mock", action="store_true", help="dry run through the same code path on TRAIN fold-2")
    ap.add_argument("--sealed-dir", type=Path, default=DEFAULT_SEALED)
    a = ap.parse_args()
    ctx = Ctx(a.mock, a.sealed_dir)
    if not SELECTED.exists():
        refuse(f"{SELECTED.relative_to(ROOT)} does not exist")
    sel = yaml.safe_load(SELECTED.read_text())
    lines = sel.get("lines", list(bm.LINES))
    check_guards(ctx, sel, lines)
    if a.mock:                                                  # fresh mock seal + output
        for d in (MOCK_SEALED, ctx.out):
            if d.exists():
                for p in d.rglob("*"):
                    if p.is_file():
                        p.chmod(0o600)
                shutil.rmtree(d)
        bm.build_mock_sealed(MOCK_SEALED, lines, ROOT)
    ctx.out.mkdir(parents=True, exist_ok=True)
    lock = {"status": "STARTED", "started": now(), "mock": a.mock, "git_hash": bm.git_hash(ROOT),
            "selected_yaml_sha256": bm.sha256_file(SELECTED)}
    ctx.lock.write_text(json.dumps(lock, indent=1))
    ctx.commit("final_eval: LOCK (STARTED) before reading any TEST file")

    labels = label_paths(ctx, lines)
    summary = {**lock, "results": [], "comparisons": [], "notes": []}
    try:
        unseal_and_verify(ctx, labels)
        token = bm.final_eval_token(ctx.sealed / "data" / "processed" / "bench", mock=a.mock)
        run_all(ctx, sel, lines, token, summary)
        lock.update(status="COMPLETED", finished=now())
    except BaseException as e:
        lock.update(status="FAILED", finished=now(), error=repr(e))
        raise
    finally:
        set_readable([p for p in labels if p.exists()], False)
        summary.update(status=lock["status"], finished=lock.get("finished"))
        (ctx.out / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
        ctx.lock.write_text(json.dumps(lock, indent=1))
        ctx.commit(f"final_eval: TEST evaluation {lock['status']}; results/final/LOCK")
    print(f"final evaluation {lock['status']} ({'mock, ' if a.mock else ''}output in {ctx.out})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
