"""Real-data invariants of the frozen benchmark tables and split (DESIGN §6, §8.6, §8.7).
Skipped when the tables have not been built. Never opens TEST label files."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from depmaptr import benchmark as bm

BD = bm.bench_dir()
SPLIT = bm.ROOT / "results" / "split"
pytestmark = pytest.mark.skipif(not (BD / "manifest.json").exists(), reason="benchmark tables not built")

# pinned expected ranges (DESIGN §8.7); built values: HCT116 5,623 x 10,209; HEK293T 2,688 x 11,000
RANGES = {"hct116": {"regs": (5000, 6500), "tgts": (9500, 11000), "pos_train": (60000, 110000)},
          "hek293t": {"regs": (2200, 3200), "tgts": (10000, 12000), "pos_train": (25000, 60000)}}


@pytest.fixture(scope="module")
def flags():
    return bm.normalise_flags(pd.read_csv(bm.ROOT / "data/processed/frozen_flags.tsv", sep="\t", dtype={"chrom": str}))


@pytest.fixture(scope="module", params=bm.LINES)
def line_tables(request):
    line = request.param
    regs = pd.read_csv(BD / f"universe_regulators_{line}.tsv", sep="\t", dtype={"chrom": str})
    tgts = pd.read_csv(BD / f"universe_targets_{line}.tsv", sep="\t", dtype={"chrom": str})
    pos = pd.read_parquet(BD / f"positives_{line}_TRAIN.parquet")
    return line, regs, tgts, pos


def test_known_arms(flags):
    sym = flags.reset_index().drop_duplicates("symbol").set_index("symbol")
    assert sym.loc["ERBB2", "arm_id"] == "17q" and sym.loc["MYC", "arm_id"] == "8q"
    assert sym.loc["CDKN2A", "arm_id"] == "9p"


def test_universe_definition(line_tables, flags):
    line, regs, tgts, pos = line_tables
    r = RANGES[line]
    assert r["regs"][0] <= len(regs) <= r["regs"][1] and r["tgts"][0] <= len(tgts) <= r["tgts"][1]
    assert regs["regulator"].is_unique and tgts["target"].is_unique and regs["pie_id"].is_unique
    f = flags.loc[regs["regulator"]]
    assert f["variable_cn"].all() and not (f["is_par"] | f["is_y"] | f["is_mt"]).any()
    assert set(regs["chrom"]) <= set(bm.CHROMS) and set(tgts["chrom"]) <= set(bm.CHROMS)
    g = flags.loc[tgts["target"]]
    assert g["measured"].all() and not (g["is_par"] | g["is_y"] | g["is_mt"]).any()
    assert (regs["chrom"].values == f["chrom"].values).all() and (regs["arm"].values == f["arm_id"].values).all()
    assert set(regs["cn_tercile"]) == {1, 2, 3}
    assert set(tgts["expr_decile"]) == set(range(1, 11))
    assert tgts["bin"].str.split(":").str[0].eq(tgts["chrom"]).all()
    # TEST label-derived columns are not in the open table
    te = regs[regs["split"] == "TEST"]
    assert te["own_class"].isna().all() and te["degree_decile"].isna().all()
    tr = regs[regs["split"] == "TRAIN"]
    assert tr["own_class"].isin(["i", "ii", "iii"]).all() and tr["fold"].isin([1, 2]).all()


def test_positives_trans_and_in_universe(line_tables):
    line, regs, tgts, pos = line_tables
    rc = dict(zip(regs["regulator"], regs["chrom"]))
    tc = dict(zip(tgts["target"], tgts["chrom"]))
    assert pos["regulator"].isin(rc).all() and pos["target"].isin(tc).all()
    assert (pos["regulator"].map(rc) != pos["target"].map(tc)).all()          # chrom(A) != chrom(B)
    assert (pos["regulator"] != pos["target"]).all()
    assert not pos.duplicated(["regulator", "target"]).any()
    assert (pos["p_adj"] < 0.05).all() and pos["pos_primary"].all()
    assert (pos.loc[pos["pos_p01"], "p_adj"] < 0.01).all() and (pos.loc[~pos["pos_p01"], "p_adj"] >= 0.01).all()
    assert (pos["ueffect"].abs().ge(0.2) & pos["pos_primary"] == pos["pos_ueffect"]).all()
    assert set(pos["lfc_sign"].unique()) <= {-1, 0, 1}
    # TRAIN positives only from TRAIN regulators
    split = dict(zip(regs["regulator"], regs["split"]))
    assert (pos["regulator"].map(split) == "TRAIN").all()
    lo, hi = RANGES[line]["pos_train"]
    assert lo <= len(pos) <= hi


def test_manifest_invariants():
    man = json.loads((BD / "manifest.json").read_text())
    for line in bm.LINES:
        assert man["lines"][line]["median_on_target_lfc_train"] < 0               # CRISPRi knocks down its target


def test_split_reproducible_and_triplets():
    units = pd.read_csv(SPLIT / "arm_units.tsv", sep="\t")
    again = bm.make_split(units[["arm", "n_hct116", "n_hek293t"]])
    assert (again[["arm", "split", "fold"]].values == units[["arm", "split", "fold"]].values).all()
    assert all((g["split"] == "TEST").sum() == 1 for _, g in units.groupby("triplet") if len(g) == 3)
    assert set(units["arm"]) == {a for a in units["arm"] if a[:-1] in bm.CHROMS}
    rs = pd.read_csv(SPLIT / "regulator_split.tsv", sep="\t", dtype={"chrom": str})
    assert rs["gene"].is_unique and rs["split"].isin(["TRAIN", "TEST"]).all()
    assert (rs["arm"].map(dict(zip(units["arm"], units["split"]))) == rs["split"]).all()
    frac = (units["split"] == "TEST").mean()
    assert 0.30 <= frac <= 0.36
    folds = pd.read_csv(SPLIT / "train_folds.tsv", sep="\t")
    assert list(folds["fold"]) == [1 + i % 2 for i in range(len(folds))]
    for line in bm.LINES:
        regs = pd.read_csv(BD / f"universe_regulators_{line}.tsv", sep="\t")
        m = regs.merge(rs, left_on="regulator", right_on="gene")
        assert len(m) == len(regs) and (m["split_x"] == m["split_y"]).all()


def test_evaluate_refuses_test_on_real_tree(tmp_path):
    with pytest.raises(bm.SealedSplitError):
        bm.evaluate(tmp_path, "hct116", "TEST")
    with pytest.raises(bm.SealedSplitError):
        bm.load_bench("hct116", "TEST")


def test_random_scores_on_real_train(tmp_path):
    """B0 smoke: random scores on the real TRAIN universe give OR ~ 1 and the pinned pair count."""
    b = bm.load_bench("hct116", "TRAIN")
    rng = np.random.default_rng(0)
    S = rng.standard_normal((len(b.regs), len(b.tgts)), dtype=np.float32)
    np.savez(tmp_path / "scores.npz", regulators=b.regs["regulator"].values.astype(str),
             targets=b.tgts["target"].values.astype(str), t=S, score=S,
             s_cis=np.ones(len(b.regs), np.int8), cis_t=np.zeros(len(b.regs), np.float32))
    # computed directly (not via evaluate(), which always logs) so that pytest runs leave no log rows
    ev, miss = bm.build_evaluation(bm.load_scores(tmp_path), b)
    m = bm.compute_metrics(ev, full=False)
    m.update(miss)
    expect = sum(len(b.tgts) - (b.tgts["chrom"] == c).sum() for c in b.regs["chrom"])
    assert m["n_pairs"] == expect
    assert 0.6 < m["primary_or"] < 1.5 and 0.85 < m["top0.01"]["or"] < 1.15
    assert m["n_regulators_missing_in_scores"] == 0
