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

selected.yaml (run ids refer to data/results/scores/<run_id>/; see results/final/selected.DRAFT.yaml):
  selected: RUN                          # the selected method
  ablations: [RUN, ...]                  # <= 3 (DESIGN §8.4)
  reported: [RUN, ...]                   # optional reported variants (e.g. the D7 context layer)
  labels: {RUN: "display label"}         # optional
  baselines: [B0, ..., B6]
  controls: [NC1, C2_g2, PC2, NC2]       # score runs, plus the built-ins PC2 and NC2
  nc2_run: RUN                           # run evaluated for NC2 (must be without G2; default: selected)
  c1: [RUN, ...]                         # D6 C1 pair-shuffle controls (default: [selected])
  sensitivity: [RUN, ...]                # e.g. selected + __S1 / __noHCT116 / ...
  subsets: [paralog_exclusion, class_i, hub_excluded, ...]   # applied to each run in subset_runs
  subset_runs: [RUN, ...]                # default: [selected]
  lines: [hct116, hek293t]               # default both
  n_perm: 1000
The older form of `subsets` entries, {name, run, subset}, is still accepted. Outputs: summary.json,
summary.tsv and summary.md (one row per role x run x line x subset x metric set; depmaptr.report).
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
from depmaptr import controls, report

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
    for sub in subset_specs(sel):
        if sub["subset"] not in KNOWN_SUBSETS:
            refuse(f"unknown subset in selected.yaml: {sub}")
    unknown = [c for c in sel.get("controls") or [] if c in BUILTIN_CONTROLS and c not in ("PC2", "NC2")]
    if unknown:
        refuse(f"unknown built-in control {unknown}")
    for r in score_runs(sel) + [x["run"] for x in subset_specs(sel)] + list(sel.get("c1") or []) + \
            [sel.get("nc2_run", sel["selected"])]:
        bm.resolve_run_dir(r)                                   # fail before the LOCK is written


BUILTIN_CONTROLS = {"PC2", "NC2"}
REPORT_ONLY_SUBSETS = {"hub_excluded"}                          # already inside every run's metrics
KNOWN_SUBSETS = set(bm.SUBSETS) | set(bm.PAIR_SUBSETS) | REPORT_ONLY_SUBSETS


def roles(sel: dict) -> list:
    """(role, run) for every score run, in evaluation order, without duplicates."""
    out, seen = [], set()
    groups = [("selected", [sel["selected"]]), ("ablation", sel.get("ablations")), ("reported", sel.get("reported")),
              ("baseline", sel.get("baselines")),
              ("control", [c for c in sel.get("controls") or [] if c not in BUILTIN_CONTROLS]),
              ("sensitivity", sel.get("sensitivity"))]
    for role, runs in groups:
        for r in runs or []:
            if r not in seen:
                seen.add(r)
                out.append((role, r))
    return out


def score_runs(sel: dict) -> list:
    return [r for _, r in roles(sel)]


all_runs = score_runs          # backwards-compatible name


def subset_specs(sel: dict) -> list:
    """Normalised subset evaluations: dicts with name, run, subset."""
    out = []
    runs = sel.get("subset_runs") or [sel["selected"]]
    for s in sel.get("subsets") or []:
        if isinstance(s, dict):
            out.append({"name": s["name"], "run": s["run"], "subset": s.get("subset", str(s["name"]).lower())})
        else:
            name = str(s).replace("-", "_")
            out += [{"name": name, "run": r, "subset": name} for r in runs]
    if "NC2" in (sel.get("controls") or []):
        out.append({"name": "NC2", "run": sel.get("nc2_run", sel["selected"]), "subset": "nc2"})
    out.append({"name": "no_split_sibling", "run": sel["selected"], "subset": "no_split_sibling"})
    return out


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
    labels = sel.get("labels") or {}
    rec = summary.setdefault("_records", [])

    def keep(role, run, line, out, subset=None):
        rec.append({"role": role, "run": run, "label": labels.get(run, ""), "line": line, "subset": subset,
                    "out": out})
        summary["results"].append({"role": role, "run": run, "subset": subset, "line": line,
                                   "split": out["split"], "primary_or": out["metrics"]["primary_or"],
                                   "perm_p": out["metrics"].get("perm_regulator", {}).get("perm_p"),
                                   "jackknife_ci95_log_or": out["metrics"]["jackknife_log_or"].get("ci95"),
                                   controls.D6_CMH_LABEL: out["metrics"].get(
                                       "cmh_regulator_x_context_exprdecile_D6", {}).get("cmh_or")})

    for line in lines:
        for role, run in roles(sel):
            keep(role, run, line, bm.evaluate(run, line, "TEST", n_perm=n_perm, _token=token))
        regs = bm.load_bench(line, "TEST", token=token).regs
        for s in subset_specs(sel):
            if s["subset"] in REPORT_ONLY_SUBSETS:
                summary["notes"].append(f"{s['name']}: reported from each run's metrics (OR_hub_excluded)")
                continue
            if s["subset"] in bm.PAIR_SUBSETS:
                out = bm.evaluate(s["run"], line, "TEST", n_perm=n_perm, subset_name=s["name"], _token=token,
                                  exclude_pairs=bm.PAIR_SUBSETS[s["subset"]](ROOT))
            else:
                mask = bm.subset_mask(regs, s["subset"], ROOT)
                out = bm.evaluate(s["run"], line, "TEST", n_perm=n_perm, reg_subset=mask, subset_name=s["name"],
                                  _token=token)
            keep("control" if s["name"] == "NC2" else "subset", s["run"], line, out, s["name"])
        for _, other in roles(sel)[1:]:
            summary["comparisons"].append(bm.compare(sel["selected"], other, line, "TEST", _token=token))
        for run in sel.get("c1") or [sel["selected"]]:               # D6 matched control C1
            keep("control", f"{run}_C1", line, controls.evaluate_c1(run, line, "TEST", n_perm=n_perm, _token=token))
            summary.setdefault("c1_comparisons", []).append(controls.compare_c1(run, line, "TEST", _token=token))
    wants_pc2 = "PC2" in (sel.get("controls") or ["PC2"])
    if wants_pc2 and {"hct116", "hek293t"} <= set(lines):
        keep("control", "PC2", "hct116", bm.evaluate_pc2("TEST", ROOT, n_perm=n_perm, _token=token))
    elif wants_pc2:
        summary["notes"].append("PC2 skipped: needs both lines")


def write_tables(ctx: Ctx, summary: dict) -> None:
    df = report.summary_frame(summary.pop("_records", []), summary.get("comparisons", []),
                              summary.get("c1_comparisons", []))
    df.to_csv(ctx.out / "summary.tsv", sep="\t", index=False)
    title = "Final TEST evaluation" + (" (MOCK: TRAIN fold-2 as pseudo-TEST)" if ctx.mock else "")
    (ctx.out / "summary.md").write_text(report.to_markdown(df, title))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mock", action="store_true", help="dry run through the same code path on TRAIN fold-2")
    ap.add_argument("--sealed-dir", type=Path, default=DEFAULT_SEALED)
    ap.add_argument("--selected", type=Path, help="(--mock only) another selection file, e.g. selected.DRAFT.yaml")
    a = ap.parse_args()
    ctx = Ctx(a.mock, a.sealed_dir)
    if a.selected is not None and not a.mock:
        refuse("--selected is only allowed with --mock; the real run uses the committed selected.yaml")
    sel_path = a.selected if a.selected is not None else SELECTED
    if not sel_path.exists():
        refuse(f"{sel_path} does not exist")
    sel = yaml.safe_load(sel_path.read_text())
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
            "selected_yaml_sha256": bm.sha256_file(sel_path), "selected_yaml": str(sel_path)}
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
        try:
            write_tables(ctx, summary)
        except Exception as e:                                   # tables must never block the LOCK update
            summary["notes"].append(f"summary tables failed: {e!r}")
        summary.pop("_records", None)
        (ctx.out / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
        ctx.lock.write_text(json.dumps(lock, indent=1))
        ctx.commit(f"final_eval: TEST evaluation {lock['status']}; results/final/LOCK")
    print(f"final evaluation {lock['status']} ({'mock, ' if a.mock else ''}output in {ctx.out})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
