#!/usr/bin/env python
"""Re-record the sha256 of changed or added in-tree frozen files in results/FROZEN.sha256 after a decision
note (DESIGN §8.5). RUN BY THE PI ONLY. Sealed entries are never touched.

Existing entries are replaced in place; new files are appended to the frozen section. A comment line
records the date, the note and the files.

Usage: .venv/bin/python scripts/benchmark/record_frozen_change.py --note "D6 controls" FILE [FILE ...]
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

from depmaptr import benchmark as bm

FROZEN = bm.ROOT / "results" / "FROZEN.sha256"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--note", required=True, help="decision id / reason (must exist in RESEARCH_TREE.md)")
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    tree = bm.ROOT / "RESEARCH_TREE.md"
    if not tree.exists() or a.note.split()[0] not in tree.read_text():
        raise SystemExit(f"decision {a.note.split()[0]!r} not found in RESEARCH_TREE.md")
    lines = FROZEN.read_text().splitlines()
    _, sealed = bm.read_frozen(FROZEN)
    sealed_paths = {q for _, q in sealed}
    cut = next(i for i, ln in enumerate(lines) if ln.startswith("# sealed"))
    head, tail = lines[:cut], lines[cut:]
    for f in a.files:
        rel = str(Path(f).resolve().relative_to(Path(bm.ROOT).resolve()))
        if rel in sealed_paths:
            raise SystemExit(f"{rel} is a sealed entry; refusing")
        entry = f"{bm.sha256_file(bm.ROOT / rel)}  {rel}"
        idx = [i for i, ln in enumerate(head) if not ln.startswith("#") and ln.split(None, 1)[-1].strip() == rel]
        if idx:
            head[idx[0]] = entry
        else:
            head.append(entry)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    head.append(f"# re-recorded {stamp} ({a.note}): {' '.join(a.files)}")
    FROZEN.write_text("\n".join(head + tail) + "\n")
    print(f"updated {len(a.files)} entries; drift now: {bm.frozen_drift(FROZEN)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
