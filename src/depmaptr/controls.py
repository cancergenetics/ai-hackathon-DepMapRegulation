"""Evaluation-only matched controls added by PI decision D6 (reviews/ideas/M5a_ideas.md §2-§3).
No frozen metric is changed; these are additional diagnostics.

C1 (pair shuffle): for each regulator A, A's score row is permuted across targets within 50 bins of the
DepMap context-expression percentile (HCT116 = ACH-000971 for the HCT116 line, the primary; HEK293 =
ACH-001085 for the HEK293T line). Percentiles are taken over ALL targets of the score file (all DepMap
measured genes), never within the PIE universe. Same-chromosome entries are shuffled only among themselves,
so masked values never move into trans pairs. C1 keeps the regulator allocation and the target
detectability, and destroys only the pair evidence.

Smoke test (HCT116 TRAIN, approximate mode): B4e, a pure target-expression score, goes from OR 7.09 to 4.08
under C1 (Δlog OR +0.55), because the top expression bin still has a steep detectability gradient. Exact mode
does not have this bias. Fused runs must therefore supply pair_score and target_weight.

If scores.npz holds the optional `pair_score` float32[R,T] and `target_weight` float32[T] (docs/INTERFACES.md),
C1 shuffles `pair_score` and recomputes score = shuffled pair_score * target_weight (exact for a fusion of
the form S = pair x target-weight, with any per-regulator factor already inside pair_score). Otherwise it
shuffles `score` itself. That is exact when the target weight is constant within a bin; for a
percentile-based weight ρ_T^γ, the within-bin variation (bin width 2 percentile points) is the
approximation error.

The permutation order comes from a seeded pair hash: SeedSequence(20261027).spawn(8)[C1_SPAWN] seeds
blake2b/splitmix64 keys per (regulator, target). It does not depend on row order or on which regulators are
evaluated.

D6 diagnostic CMH: regulator x context-expression decile (deciles of DepMap context expression among the
line's universe targets). It is reported next to the frozen regulator x S2-mean-expression CMH, which is
unchanged.
"""
from __future__ import annotations

import datetime as _dt
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from depmaptr import benchmark as bm

C1_SPAWN = 4                      # benchmark.SPAWN uses 0-3
C1_BINS = 50
D6_CMH_LABEL = "CMH_OR(reg x HCT116 expr) [D6 diagnostic]"
CONTEXT_COL = {"hct116": "expr_hct116", "hek293t": "expr_hek293"}


def context_expression(line: str, targets, root: Path = bm.ROOT) -> np.ndarray:
    """DepMap log2(TPM+1) of `targets` in the line's context model (NaN where absent)."""
    p = Path(root) / "data" / "processed" / "depmap_expr_matched.parquet"
    row = pd.read_parquet(p, filters=[("ModelID", "==", bm.CONTEXT_MODEL[line])])
    if len(row) != 1:
        raise ValueError(f"context model {bm.CONTEXT_MODEL[line]} not found once in {p}")
    s = row.drop(columns="ModelID").iloc[0]
    return s.reindex(pd.Index(targets)).to_numpy(np.float64)


def c1_bins(expr: np.ndarray, n_bins: int = C1_BINS) -> np.ndarray:
    """0-based percentile bins over all score-file targets (average ranks; NaN expression -> extra bin)."""
    b = np.full(len(expr), n_bins, dtype=np.int64)
    ok = np.isfinite(expr)
    b[ok] = bm.quantile_bins(expr[ok], n_bins) - 1
    return b


def _shuffle_rows(M: np.ndarray, reg_keys, tgt_keys, reg_chrom, tgt_chrom, bins: np.ndarray,
                  seed: int, chunk: int = 256) -> np.ndarray:
    """Permute each row of M within (bin x trans/cis) groups using seeded pair-hash keys."""
    ss = np.random.SeedSequence(seed).spawn(bm.N_SPAWN)[C1_SPAWN]
    hseed = int(ss.generate_state(1)[0])
    hr = bm.key_hash(reg_keys, hseed)
    ht = bm.key_hash(tgt_keys, hseed)
    out = np.empty_like(M)
    nb = int(bins.max()) + 1
    tgt_chrom = np.asarray(tgt_chrom).astype(str)
    reg_chrom = np.asarray(reg_chrom).astype(str)
    for c in np.unique(reg_chrom):
        rows = np.flatnonzero(reg_chrom == c)
        grp = bins + np.where(tgt_chrom == c, nb, 0)             # cis targets form separate groups
        canon = np.argsort(grp, kind="stable")                  # slot k of the bin-sorted order
        for i in range(0, len(rows), chunk):
            r = rows[i:i + chunk]
            keys = (bm.pair_hash(hr[r][:, None], ht[None, :]) >> np.uint64(11)).astype(np.float64) * 2.0 ** -53
            order = np.argsort(grp[None, :] + keys, axis=1)     # within-group random order
            blk = M[r]
            out[np.ix_(r, canon)] = np.take_along_axis(blk, order, axis=1)
    return out


def c1_scores(sc: dict, line: str, regulators=None, root: Path = bm.ROOT, seed: int = bm.SEED,
              chrom: dict | None = None) -> dict:
    """C1 version of a score dict (from benchmark.load_scores), restricted to `regulators` if given.
    `chrom`: gene key -> chromosome (default: frozen flags)."""
    R = np.asarray(sc["regulators"]).astype(str)
    T = np.asarray(sc["targets"]).astype(str)
    if regulators is not None:
        keep = np.flatnonzero(pd.Index(R).isin(list(regulators)))
    else:
        keep = np.arange(len(R))
    if chrom is None:
        f = bm.normalise_flags(pd.read_csv(Path(root) / "data/processed/frozen_flags.tsv", sep="\t",
                                           dtype={"chrom": str}))
        chrom = f["chrom"].to_dict()
    rc = np.array([chrom.get(g, "?") for g in R[keep]])
    tc = np.array([chrom.get(g, "?") for g in T])
    bins = c1_bins(context_expression(line, T, root))
    exact = "pair_score" in sc and "target_weight" in sc
    src = np.asarray(sc["pair_score"] if exact else sc["score"])[keep]
    sh = _shuffle_rows(src, R[keep], T, rc, tc, bins, seed)
    score = (sh * np.asarray(sc["target_weight"], np.float32)[None, :]).astype(np.float32) if exact else sh
    meta = dict(sc.get("meta", {}))
    base = meta.get("run_id", "scores")
    meta.update(run_id=f"{base}_C1", variant=f"{meta.get('variant', base)}+C1",
                params={**meta.get("params", {}), "control": "C1", "c1_bins": C1_BINS, "c1_line": line,
                        "c1_mode": "pair_score x target_weight" if exact else "score (approximation)"})
    return {"regulators": R[keep], "targets": T, "score": score, "s_cis": np.asarray(sc["s_cis"])[keep],
            "meta": meta, "c1_exact": exact}


def d6_cmh(ev: "bm.Evaluation", sel, pl, line: str) -> dict:
    """Regulator x context-expression decile CMH (deciles among the evaluated universe targets)."""
    col = CONTEXT_COL[line]
    te = bm.quantile_bins(ev.tgts[col].astype(float).fillna(-1.0).values, 10) - 1
    saved = ev.te
    try:
        ev.te = te
        out = bm.mantel_haenszel(*ev.cmh_tables(sel, pl, "regulator"))
    finally:
        ev.te = saved
    out["label"] = D6_CMH_LABEL if line == "hct116" else "CMH_OR(reg x HEK293 expr) [D6 diagnostic]"
    return out


def _paired_jackknife(ev_a, ev_b, sel_a, sel_b, pl) -> dict:
    """Delete-one-arm jackknife (t_{G-1}) of Δlog OR and of ΔOR = OR_a − OR_b; selections fixed."""
    out = {}
    ca, cb = ev_a.arm_counts(sel_a, pl), ev_b.arm_counts(sel_b, pl)
    G = ev_a.G

    def ors(c, drop=None):
        a, s, p, n = c
        if drop is None:
            return bm.odds_ratio(*bm.table_from_counts(a.sum(), s.sum(), p.sum(), n.sum()))
        return bm.odds_ratio(*bm.table_from_counts(a.sum() - a, s.sum() - s, p.sum() - p, n.sum() - n))

    full_a, full_b = float(ors(ca)), float(ors(cb))
    del_a, del_b = ors(ca, True), ors(cb, True)
    tq = stats.t.ppf(0.975, G - 1)
    for name, est, d in (("delta_log_or", math.log(full_a) - math.log(full_b), np.log(del_a) - np.log(del_b)),
                         ("delta_or", full_a - full_b, del_a - del_b)):
        se = math.sqrt((G - 1) / G * np.sum((d - d.mean()) ** 2))
        out[name] = {"estimate": est, "se": se, "ci95": [est - tq * se, est + tq * se],
                     "lb90_one_sided": est - stats.t.ppf(0.95, G - 1) * se}
    out.update(or_run=full_a, or_c1=full_b, G=G)
    return out


def compare_c1(run, line: str = bm.PRIMARY_LINE, split: str = "TRAIN", *, n_boot: int = bm.N_BOOT,
               root: Path = bm.ROOT, debug: bool = False, _token=None) -> dict:
    """Pair-specific gain (D6): OR(run) − OR(run under C1), with a paired arm jackknife and a paired two-way
    bootstrap of Δlog OR. Both scores share the run's top-K convention. Always logged (kind control:C1)."""
    t0 = time.time()
    root = Path(root)
    bench = bm.load_bench(line, split, root, token=_token)
    sc, checks = bm._resolve_scores(run, root)
    sc1 = c1_scores(sc, line, bench.regs["regulator"], root)
    ev, _ = bm.build_evaluation(sc, bench)
    ev1, _ = bm.build_evaluation(sc1, bench)
    pl = ev.pos_lin("pos_primary")
    sel, sel1 = ev.topk(bm.PRIMARY_K), ev1.topk(bm.PRIMARY_K)
    res = {"run": sc["meta"].get("run_id"), "control": "C1", "line": line, "split": bench.split,
           "c1_mode": sc1["meta"]["params"]["c1_mode"], **_paired_jackknife(ev, ev1, sel, sel1, pl)}
    if not sc1["c1_exact"]:
        res["c1_warning"] = ("approximate C1: `score` shuffled directly. This is biased upward for scores whose "
                             "target factor varies within bins (smoke test: B4e, a target-only score, shows a "
                             "spurious Δlog OR of +0.55). Supply pair_score and target_weight for fused runs.")
    W = ev.boot_weights(n_boot)
    bd = ev.bootstrap_log_or(sel, pl, W) - ev1.bootstrap_log_or(sel1, pl, W)
    res["boot_delta_log_or_lb90"] = float(np.quantile(bd, 0.05))
    res["boot_delta_log_or_ci95"] = [float(np.quantile(bd, .025)), float(np.quantile(bd, .975))]
    res["d6_cmh_run"] = d6_cmh(ev, sel, pl, line)["cmh_or"]
    res["d6_cmh_c1"] = d6_cmh(ev1, sel1, pl, line)["cmh_or"]
    res["wall_s"] = time.time() - t0
    res = bm._jsonable(res)
    bm._append_log({"timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
                    "run_id": f"{res['run']}|C1", "variant": str(sc["meta"].get("variant", "")),
                    "kind": "control:C1", "line": line, "split": bench.split,
                    "params_hash": bm._params_hash(sc["meta"], {"control": "C1", "bins": C1_BINS}),
                    "git_hash": bm.git_hash(root),
                    "data_checksums": json.dumps({Path(k).name: v[:16] for k, v in {**bench.files, **checks}.items()}),
                    "seed": bm.SEED, "spawn": json.dumps({"c1": C1_SPAWN, **bm.SPAWN}), "n_perm": 0,
                    "n_boot": n_boot, "metrics": json.dumps(res), "wall_s": f"{res['wall_s']:.2f}"}, root, debug)
    return res


def evaluate_c1(run, line: str = bm.PRIMARY_LINE, split: str = "TRAIN", *, root: Path = bm.ROOT,
                debug: bool = False, _token=None, **kw) -> dict:
    """Full evaluate() of the C1-shuffled score (logged as run_id <run>_C1)."""
    bench = bm.load_bench(line, split, root, token=_token)
    sc, _ = bm._resolve_scores(run, root)
    return bm.evaluate(c1_scores(sc, line, bench.regs["regulator"], root), line, split, root=root,
                       debug=debug, _token=_token, **kw)
