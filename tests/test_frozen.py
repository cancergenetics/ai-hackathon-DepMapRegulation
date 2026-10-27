"""DESIGN §8.5: fail on any drift of the frozen files recorded in results/FROZEN.sha256, and check that the
sealed files are no longer in the project tree. Skipped until the freeze."""
from __future__ import annotations

import pytest

from depmaptr import benchmark as bm

FROZEN = bm.ROOT / "results" / "FROZEN.sha256"


def entries():
    frozen, sealed, section = [], [], "frozen"
    for ln in FROZEN.read_text().splitlines():
        if ln.startswith("# sealed"):
            section = "sealed"
        elif ln.startswith("# FROZEN") or ln.startswith("# frozen"):
            section = "frozen"
        if not ln.strip() or ln.startswith("#"):
            continue
        h, p = ln.split(None, 1)
        (sealed if section == "sealed" else frozen).append((h, p.strip()))
    return frozen, sealed


@pytest.mark.skipif(not FROZEN.exists(), reason="not frozen yet")
def test_no_drift_in_frozen_files():
    frozen, _ = entries()
    assert frozen
    drift = [p for h, p in frozen if not (bm.ROOT / p).exists() or bm.sha256_file(bm.ROOT / p) != h]
    assert not drift, f"frozen files changed (needs a decision note, DESIGN §8.5): {drift}"


@pytest.mark.skipif(not FROZEN.exists(), reason="not frozen yet")
def test_sealed_files_not_in_tree():
    _, sealed = entries()
    present = [p for _, p in sealed if (bm.ROOT / p).exists()]
    assert not present, f"sealed files are back in the tree: {present}"
