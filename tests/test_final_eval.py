"""Integration test of the TEST-sealing workflow on a throwaway toy project (never on the real data):
final_eval.py --mock (same sealed-file code path on TRAIN fold-2) -> seal_test.sh -> evaluate() still works on
TRAIN, TEST unreadable -> drift / uncommitted guards -> final_eval.py commits LOCK (STARTED) before reading TEST,
then COMPLETED -> any rerun refuses; a crash after the LOCK still consumes the run."""
from __future__ import annotations

import json

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
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
    (root / "results" / "split").mkdir(parents=True, exist_ok=True)
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
    (fin / "selected.yaml").write_text(yaml.safe_dump({
        "selected": "good", "ablations": ["random"], "reported": ["good"], "labels": {"random": "toy ablation"},
        "baselines": ["random"], "controls": ["PC2", "NC2"], "nc2_run": "random", "c1": ["good"],
        "subsets": ["paralog_exclusion", "class_i", "hub-excluded"], "lines": ["hct116", "hek293t"],
        "n_perm": 19}))
    regs = pd.read_csv(root / "data/processed/bench/universe_regulators_hct116.tsv", sep="\t")
    tg = pd.read_csv(root / "data/processed/bench/universe_targets_hct116.tsv", sep="\t")
    pd.DataFrame({"gene_a": list(regs.regulator[:6]), "gene_b": list(tg.target[:6])}).to_csv(
        root / "data/processed/paralog_pairs.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": regs.regulator, "reg_class": ["TF", "other"] * (len(regs) // 2) + ["other"] * (len(regs) % 2)}
                 ).to_csv(root / "data/processed/regulator_classes.tsv", sep="\t", index=False)
    sealed = tmp_path / "sealed"
    run(["git", "init", "-q"], root, sealed, check=True)
    run(["git", "add", "results/final/selected.yaml"], root, sealed, check=True)
    run(["git", "commit", "-qm", "selected"], root, sealed, check=True)
    return root, sealed


def final_eval(root, sealed, *args):
    return run([sys.executable, str(REPO / "scripts" / "final_eval.py"), "--sealed-dir", str(sealed), *args],
               root, sealed)


def lock(root, mock=False):
    return json.loads((root / "results" / ("final_mock" if mock else "final") / "LOCK").read_text())


def test_seal_and_final_eval_workflow(toy_project):
    root, sealed = toy_project
    # a real run before the seal is refused
    r = final_eval(root, sealed)
    assert r.returncode != 0 and "FROZEN.sha256 missing" in r.stderr
    # dry run (DESIGN §6) through the sealed-file code path: mock seal of TRAIN fold-2 labels
    r = final_eval(root, sealed, "--mock")
    assert r.returncode == 0, r.stderr
    assert lock(root, mock=True)["status"] == "COMPLETED"
    summ = json.loads((root / "results" / "final_mock" / "summary.json").read_text())
    assert {x["run"] for x in summ["results"]} == {"good", "random", "PC2", "good_C1"}
    assert len(summ["c1_comparisons"]) == 2 and all(c["split"] == "PSEUDO_TEST" for c in summ["c1_comparisons"])
    from depmaptr.controls import D6_CMH_LABEL
    assert all(x[D6_CMH_LABEL] is not None for x in summ["results"] if x["run"] == "good" and not x["subset"])
    assert {x.get("subset") for x in summ["results"]} >= {"NC2", "no_split_sibling", "paralog_exclusion", "class_i"}
    assert any("hub_excluded" in n for n in summ["notes"])
    tab = pd.read_csv(root / "results" / "final_mock" / "summary.tsv", sep="\t")
    assert {"pooled", "cmh", "sign", "C1", "paired", "runtime", "class:TF"} <= set(tab.metric_set)
    assert set(tab.line) == {"hct116", "hek293t"} and (tab.split.dropna() == "PSEUDO_TEST").all()
    assert (tab[(tab.role == "ablation") & (tab.metric_set == "pooled")].label == "toy ablation").all()
    md = (root / "results" / "final_mock" / "summary.md").read_text()
    assert "## pooled" in md and "## C1" in md and "MOCK" in md
    assert all(x["split"] == "PSEUDO_TEST" for x in summ["results"])
    for q in (root / "data" / "interim" / "mock_sealed").rglob("*_TEST.*"):
        with pytest.raises(PermissionError):                # chmod round trip restored the seal
            open(q, "rb")
    assert not (root / "results" / "final" / "LOCK").exists()
    assert final_eval(root, sealed, "--mock").returncode == 0           # the mock can be repeated

    r = run(["bash", "scripts/benchmark/seal_test.sh"], root, sealed)
    assert r.returncode == 0, r.stderr
    for rel in ("data/processed/bench/positives_hct116_TEST.parquet",
                "data/processed/bench/regulator_labels_hct116_TEST.tsv",
                "data/processed/bench/positives_hek293t_TEST.parquet",
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
    assert run(["bash", "scripts/benchmark/seal_test.sh"], root, sealed).returncode != 0   # only once
    # TRAIN evaluation still works after the seal; TEST stays refused
    out = bm.evaluate(root / "data/results/scores/good", "hct116", "TRAIN", n_perm=9, n_boot=9, root=root)
    assert out["metrics"]["primary_or"] > 1
    with pytest.raises(bm.SealedSplitError):
        bm.evaluate(root / "data/results/scores/good", "hct116", "TEST", root=root)

    # drift in a frozen in-tree file blocks the real run (before any LOCK)
    tf = root / "data" / "processed" / "bench" / "universe_targets_hct116.tsv"
    orig = tf.read_bytes()
    tf.write_bytes(orig + b"\n")
    r = final_eval(root, sealed)
    assert r.returncode != 0 and "frozen in-tree files changed" in r.stderr
    assert not (root / "results" / "final" / "LOCK").exists()
    tf.write_bytes(orig)
    # an uncommitted change to selected.yaml blocks the real run
    sel = root / "results" / "final" / "selected.yaml"
    sel.write_text(sel.read_text() + "\n# edit\n")
    r = final_eval(root, sealed)
    assert r.returncode != 0 and "must be committed" in r.stderr
    run(["git", "checkout", "--", "results/final/selected.yaml"], root, sealed, check=True)

    # the one real run
    r = final_eval(root, sealed)
    assert r.returncode == 0, r.stderr
    assert lock(root)["status"] == "COMPLETED"
    msgs = run(["git", "log", "--format=%s"], root, sealed).stdout.splitlines()
    assert msgs[0].startswith("final_eval: TEST evaluation COMPLETED")
    assert msgs[1].startswith("final_eval: LOCK (STARTED)")              # committed before TEST was read
    first = run(["git", "show", "--name-only", "--format=", "HEAD~1"], root, sealed).stdout.split()
    assert first == ["results/final/LOCK"]
    assert (root / "results" / "final" / "summary.json").exists()
    with pytest.raises(PermissionError):
        open(sealed / "data/processed/bench/positives_hct116_TEST.parquet", "rb")
    r = final_eval(root, sealed)                                          # any rerun refuses
    assert r.returncode != 0 and "LOCK" in r.stderr


def test_final_eval_refuses_without_selected(toy_project):
    root, sealed = toy_project
    (root / "results" / "final" / "selected.yaml").unlink()
    r = final_eval(root, sealed, "--mock")
    assert r.returncode != 0 and "selected.yaml" in r.stderr


def test_crash_after_lock_consumes_the_run(toy_project):
    """A tampered sealed file is detected after the STARTED LOCK is committed: the run is FAILED, committed,
    and cannot be retried."""
    root, sealed = toy_project
    assert run(["bash", "scripts/benchmark/seal_test.sh"], root, sealed).returncode == 0
    p = sealed / "data/processed/bench/regulator_labels_hct116_TEST.tsv"
    p.chmod(0o600)
    p.write_text(p.read_text() + "X\t0\t0\n")
    p.chmod(0o000)
    r = final_eval(root, sealed)
    assert r.returncode != 0 and "sha256 mismatch" in r.stderr
    lk = lock(root)
    assert lk["status"] == "FAILED" and "sha256 mismatch" in lk["error"]
    assert run(["git", "cat-file", "-e", "HEAD:results/final/LOCK"], root, sealed).returncode == 0
    r = final_eval(root, sealed)
    assert r.returncode != 0 and "LOCK" in r.stderr


def test_record_frozen_change(toy_project):
    root, sealed = toy_project
    assert run(["bash", "scripts/benchmark/seal_test.sh"], root, sealed).returncode == 0
    tf = root / "data" / "processed" / "bench" / "universe_targets_hct116.tsv"
    tf.write_bytes(tf.read_bytes() + b"\n")
    script = [sys.executable, str(REPO / "scripts" / "benchmark" / "record_frozen_change.py")]
    rel = "data/processed/bench/universe_targets_hct116.tsv"
    r = run(script + ["--note", "D9 test", rel], root, sealed)
    assert r.returncode != 0 and "not found in RESEARCH_TREE" in r.stderr     # needs a decision note
    (root / "RESEARCH_TREE.md").write_text("- D9: test change\n")
    before = bm.read_frozen(root / "results" / "FROZEN.sha256")[1]
    r = run(script + ["--note", "D9 test", rel], root, sealed)
    assert r.returncode == 0, r.stderr
    assert bm.frozen_drift(root / "results" / "FROZEN.sha256", root) == []
    assert bm.read_frozen(root / "results" / "FROZEN.sha256")[1] == before          # sealed block untouched
    r = run(script + ["--note", "D9", "data/processed/bench/positives_hct116_TEST.parquet"], root, sealed)
    assert r.returncode != 0
