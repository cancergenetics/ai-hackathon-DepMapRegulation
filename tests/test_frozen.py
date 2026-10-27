"""DESIGN §8.5: fail on any drift of the frozen files recorded in results/FROZEN.sha256, and check that the
sealed files are no longer in the project tree. Skipped until the freeze."""
from __future__ import annotations

import pytest

from depmaptr import benchmark as bm

FROZEN = bm.ROOT / "results" / "FROZEN.sha256"


def entries():
    return bm.read_frozen(FROZEN)


@pytest.mark.skipif(not FROZEN.exists(), reason="not frozen yet")
def test_no_drift_in_frozen_files():
    frozen, _ = entries()
    assert frozen
    drift = bm.frozen_drift(FROZEN)
    assert not drift, f"frozen files changed (needs a decision note, DESIGN §8.5): {drift}"


@pytest.mark.skipif(not FROZEN.exists(), reason="not frozen yet")
def test_sealed_files_not_in_tree():
    _, sealed = entries()
    present = [p for _, p in sealed if (bm.ROOT / p).exists()]
    assert not present, f"sealed files are back in the tree: {present}"
