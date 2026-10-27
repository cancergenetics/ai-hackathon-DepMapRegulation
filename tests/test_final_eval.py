"""Integration test of the TEST-sealing workflow on a throwaway toy project (never on the real data):
seal_test.sh -> evaluate() still works on TRAIN, TEST unreadable -> final_eval.py --mock -> final_eval.py
(real) writes and commits LOCK -> any rerun refuses."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from depmaptr import benchmark as bm
from test_benchmark import write_toy_bench

REPO = Path(__file__).resolve().parents[1]


def run(cmd, root, sealed, check=False):
    env = dict(os.environ, DEPMAPTR_ROOT=str(root), DEPMAPTR_SEALED=str(sealed), GIT_AUTHOR_NAME="t",
               GIT_AUTHOR_EMAIL="t@t", GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t")
    r = subprocess.run(cmd, cwd=root, env=env, capture_output=True, text=True)
    if check and r.returncode != 0:
        raise AssertionError(r.stdout + r.stderr)
    return r


@pytest.fixture
def toy_project(tmp_path):
    root = tmp_path / "proj"
    root.mkdir()
    write_toy_bench(root)
    (root / "results" / "split").mkdir(parents=True)
    (root / "results" / "split" / "regulator_split.tsv").write_text("gene\tarm\tsplit\n")
    (root / "data" / "processed" / "bench" / "manifest.json").write_text("{}")
    de = root / "data" / "external" / "pie_x_atlas_orion" / "de"
    de.mkdir(parents=True)
    (de / "hct116.parquet").write_bytes(b"raw DE table")
    pie = root / "data" / "processed" / "pie"
    pie.mkdir(parents=True)
    for f in ("tested_hct116.parquet", "loose_hct116.parquet", "targets_hct116.tsv", "vocab_targets.tsv",
              "features_hct116.tsv"):
        (pie / f).write_text(f)
    (root / "scripts" / "benchmark").mkdir(parents=True)
    shutil.copy(REPO / "scripts" / "benchmark" / "seal_test.sh", root / "scripts" / "benchmark")
    fin = root / "results" / "final"
    fin.mkdir(parents=True)
    (fin / "selected.yaml").write_text(yaml.safe_dump({"selected": "good", "baselines": ["random"],
                                                       "lines": ["hct116"], "n_perm": 19}))
    sealed = tmp_path / "sealed"
    run(["git", "init", "-q"], root, sealed, check=True)
    run(["git", "add", "results/final/selected.yaml"], root, sealed, check=True)
    run(["git", "commit", "-qm", "selected"], root, sealed, check=True)
    return root, sealed


def final_eval(root, sealed, *args):
    return run([sys.executable, str(REPO / "scripts" / "final_eval.py"), "--sealed-dir", str(sealed), *args],
               root, sealed)


def test_seal_and_final_eval_workflow(toy_project):
    root, sealed = toy_project
    # dry run first (DESIGN §6): pseudo-TEST from TRAIN fold-2 arms, nothing unsealed, no LOCK
    r = final_eval(root, sealed, "--mock")
    assert r.returncode == 0, r.stderr
    assert (root / "results" / "final_mock" / "summary.json").exists()
    assert not (root / "results" / "final" / "LOCK").exists()

    r = run(["bash", "scripts/benchmark/seal_test.sh"], root, sealed)
    assert r.returncode == 0, r.stderr
    bench = root / "data" / "processed" / "bench"
    for rel in ("data/processed/bench/positives_hct116_TEST.parquet",
                "data/processed/bench/regulator_labels_hct116_TEST.tsv",
                "data/external/pie_x_atlas_orion/de/hct116.parquet", "data/processed/pie/tested_hct116.parquet",
                "data/processed/pie/loose_hct116.parquet", "data/processed/pie/targets_hct116.tsv"):
        assert not (root / rel).exists()
        assert (sealed / rel).exists()
        with pytest.raises(PermissionError):
            open(sealed / rel, "rb")
    for rel in ("data/processed/pie/vocab_targets.tsv", "data/processed/pie/features_hct116.tsv",
                "data/processed/bench/positives_hct116_TRAIN.parquet"):
        assert (root / rel).exists()
    frozen = (root / "results" / "FROZEN.sha256").read_text()
    assert "# sealed" in frozen and "positives_hct116_TEST.parquet" in frozen
    assert "positives_hct116_TRAIN.parquet" in frozen
    # sealing twice is refused
    assert run(["bash", "scripts/benchmark/seal_test.sh"], root, sealed).returncode != 0
    # TRAIN evaluation still works after the seal; TEST stays refused
    out = bm.evaluate(root / "data/results/scores/good", "hct116", "TRAIN", n_perm=9, n_boot=9, root=root)
    assert out["metrics"]["primary_or"] > 1
    with pytest.raises(bm.SealedSplitError):
        bm.evaluate(root / "data/results/scores/good", "hct116", "TEST", root=root)

    # an uncommitted change to selected.yaml blocks the real run
    sel = root / "results" / "final" / "selected.yaml"
    sel.write_text(sel.read_text() + "\n# edit\n")
    r = final_eval(root, sealed)
    assert r.returncode != 0 and "must be committed" in r.stderr
    run(["git", "checkout", "--", "results/final/selected.yaml"], root, sealed, check=True)

    # the one real run: verifies hashes, unseals, evaluates, re-seals, writes and commits LOCK
    r = final_eval(root, sealed)
    assert r.returncode == 0, r.stderr
    assert (root / "results" / "final" / "LOCK").exists()
    assert (root / "results" / "final" / "summary.json").exists()
    assert "LOCK" in run(["git", "show", "--stat", "HEAD"], root, sealed).stdout
    with pytest.raises(PermissionError):
        open(sealed / "data/processed/bench/positives_hct116_TEST.parquet", "rb")
    # any rerun refuses
    r = final_eval(root, sealed)
    assert r.returncode != 0 and "LOCK" in r.stderr


def test_final_eval_refuses_without_selected(toy_project):
    root, sealed = toy_project
    (root / "results" / "final" / "selected.yaml").unlink()
    r = final_eval(root, sealed, "--mock")
    assert r.returncode != 0 and "selected.yaml" in r.stderr


def test_final_eval_refuses_tampered_sealed_file(toy_project):
    root, sealed = toy_project
    assert run(["bash", "scripts/benchmark/seal_test.sh"], root, sealed).returncode == 0
    p = sealed / "data/processed/bench/regulator_labels_hct116_TEST.tsv"
    p.chmod(0o600)
    p.write_text(p.read_text() + "X\t0\t0\n")
    p.chmod(0o000)
    r = final_eval(root, sealed)
    assert r.returncode != 0 and "sha256 mismatch" in r.stderr
    assert not (root / "results" / "final" / "LOCK").exists()
