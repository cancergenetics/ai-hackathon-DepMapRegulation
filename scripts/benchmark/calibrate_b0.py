#!/usr/bin/env python
"""Calibration of the harness's inference on the real TRAIN universe under B0 (random scores; DESIGN §6:
"The permutation tests are calibrated on B0 and NC1 (>= 200 seeds: the p-values must be ~ uniform)").

For each seed, scores ~ N(0, 1) i.i.d. over the TRAIN universe (score sign random, s_cis = +1). Records the
regulator-label permutation p (own-gene x CN tercile), the degree-decile permutation p, the Fisher p, the
CMH p, the sign-permutation p, and whether the arm-jackknife 95 % CI and the two-way bootstrap 95 % CI of
log OR cover 0 (true log OR = 0 under B0).
Seeds: numpy.random.SeedSequence(20261027).spawn(100 + n)[100 + s] for s = 0..n-1 (spawn indices 100+).
TRAIN only. Output: results/benchmark_calibration/b0_<line>_TRAIN.tsv (+ a summary JSON).

Usage: .venv/bin/python scripts/benchmark/calibrate_b0.py [--n-seeds 200] [--line hct116] [--n-perm 1000]
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

from depmaptr import benchmark as bm

SPAWN_OFFSET = 100


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-seeds", type=int, default=200)
    ap.add_argument("--line", default="hct116", choices=bm.LINES)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--n-boot", type=int, default=200)
    a = ap.parse_args()
    bench = bm.load_bench(a.line, "TRAIN")
    regs, tgts, pos = bench.regs, bench.tgts, bench.pos
    children = np.random.SeedSequence(bm.SEED).spawn(SPAWN_OFFSET + a.n_seeds)
    rows = []
    t0 = time.time()
    for s in range(a.n_seeds):
        rng = np.random.default_rng(children[SPAWN_OFFSET + s])
        S = rng.standard_normal((len(regs), len(tgts)), dtype=np.float32)
        ev = bm.Evaluation(regs, tgts, pos, S, np.ones(len(regs), np.int8))
        sel = ev.topk(bm.PRIMARY_K)
        pl = ev.pos_lin()
        r = ev.pooled(sel, pl)
        st1 = (ev.regs["own_class"].astype(str) + "|" + ev.regs["cn_tercile"].astype(str)).values
        st2 = (ev.regs["own_class"].astype(str) + "|" + ev.regs["degree_decile"].astype(str)).values
        p1 = ev.permutation(sel, pl, st1, a.n_perm, rng)
        p2 = ev.permutation(sel, pl, st2, a.n_perm, rng)
        jk = ev.jackknife([sel], [pl])
        bl = ev.bootstrap_log_or(sel, pl, ev.boot_weights(a.n_boot, rng))
        lo, hi = np.quantile(bl, [0.025, 0.975])
        cm = ev.cmh(sel, pl)
        sg = ev.sign_metrics(sel, n_perm=200, rng=rng)
        rows.append({"seed_index": SPAWN_OFFSET + s, "or": r["or"], "a": r["a"],
                     "fisher_p": bm.fisher(r["a"], r["b"], r["c"], r["d"])["fisher_p_greater"],
                     "perm_p": p1["perm_p"], "perm_degree_p": p2["perm_p"],
                     "median_null_or": p1["median_null_or"], "jk_covers_0": jk["ci95"][0] <= 0 <= jk["ci95"][1],
                     "jk_se": jk["se"], "boot_covers_0": lo <= 0 <= hi, "cmh_p": cm.get("cmh_p"),
                     "sign_perm_p": sg.get("sign_perm_p")})
        if s % 20 == 0:
            print(f"seed {s} {time.time() - t0:.0f}s", file=sys.stderr)
    df = pd.DataFrame(rows)
    od = bm.ROOT / "results" / "benchmark_calibration"
    od.mkdir(parents=True, exist_ok=True)
    df.to_csv(od / f"b0_{a.line}_TRAIN.tsv", sep="\t", index=False)
    summ = {"line": a.line, "n_seeds": a.n_seeds, "n_perm": a.n_perm, "n_boot": a.n_boot,
            "wall_s": time.time() - t0}
    for c in ("perm_p", "perm_degree_p", "fisher_p", "cmh_p", "sign_perm_p"):
        v = df[c].dropna().values
        summ[c] = {"mean": float(v.mean()), "frac_le_0.05": float(np.mean(v <= 0.05)),
                   "frac_le_0.10": float(np.mean(v <= 0.10)), "ks_uniform_p": float(stats.kstest(v, "uniform").pvalue)}
    summ["jk_ci95_coverage"] = float(df["jk_covers_0"].mean())
    summ["boot_ci95_coverage"] = float(df["boot_covers_0"].mean())
    summ["median_or"] = float(df["or"].median())
    (od / f"b0_{a.line}_TRAIN_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
