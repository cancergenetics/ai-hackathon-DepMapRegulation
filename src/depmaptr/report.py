"""Compact result tables (D8): one row per (role, run, line, subset) x metric set, from evaluate() outputs.

Used by scripts/final_eval.py for results/final/summary.tsv and summary.md, and usable on TRAIN outputs.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from depmaptr.controls import D6_CMH_LABEL


def _exp(x):
    return None if x is None else float(math.exp(x))


def metric_sets(out: dict) -> dict:
    """metric set -> {column: value} for one evaluate() output."""
    m = out["metrics"]
    k = m["top0.001"]
    jk = m.get("jackknife_log_or", {})
    sets = {"pooled": {"n_pairs": m["n_pairs"], "K": k["K"], "OR@0.1%": k["or"], "OR@0.01%": m["top0.0001"]["or"],
                       "OR@1%": m["top0.01"]["or"], "fisher_p(indep.)": k["fisher_p_greater"],
                       "jk_CI95_lo": _exp((jk.get("ci95") or [None, None])[0]),
                       "jk_CI95_hi": _exp((jk.get("ci95") or [None, None])[1])}}
    if "perm_regulator" in m:
        sets["pooled"].update({"perm_p": m["perm_regulator"]["perm_p"], "perm_p_degree": m["perm_degree"]["perm_p"],
                               "OR/median_null": m["perm_regulator"]["or_over_median_null"],
                               "OR/median_null_degree": m["perm_degree"]["or_over_median_null"]})
        sets["cmh"] = {"CMH_OR(reg x S2 expr)": m["cmh_regulator_x_exprdecile"]["cmh_or"],
                       "CMH_OR(own class x S2 expr)": m["cmh_ownclass_x_exprdecile"]["cmh_or"],
                       D6_CMH_LABEL: m.get("cmh_regulator_x_context_exprdecile_D6", {}).get("cmh_or"),
                       "mean_AUROC": m["auroc"]["mean_auroc"], "top10_share": m["top10_regulator_share"],
                       "OR_hub_excluded": m["hub_excluded"]["or"]}
        sg = m["sign"]
        sets["sign"] = {"sign_concordant_OR": sg["sign_concordant_or"], "n_overlap": sg["n_overlap"],
                        "concordance": sg.get("concordance"), "c0": sg.get("c0"), "sign_OR": sg.get("sign_or"),
                        "sign_perm_p": sg.get("sign_perm_p"),
                        "reg_weighted_concordance": sg.get("regulator_weighted_concordance")}
    if "callset" in m and "or" in m["callset"]:
        sets["callset"] = {"n_calls": m["callset"]["n_calls"], "callset_OR": m["callset"]["or"],
                           "callset_fisher_p": m["callset"]["fisher_p_greater"]}
    for cls, r in (m.get("regulator_class") or {}).items():
        sets[f"class:{cls}"] = {"n_regulators": r["n_regulators"], "selected": r["selected"],
                                "positives": r["positives"], "a": r["a"], "OR@0.1%": r["or"],
                                "fisher_p(indep.)": r["fisher_p_greater"]}
    sets["runtime"] = {"eval_wall_s": out.get("wall_s"), "score_wall_s": out.get("score_wall_s"),
                       "threads": out.get("threads")}
    return sets


def summary_frame(records: list, comparisons: list = (), c1: list = ()) -> pd.DataFrame:
    """records: dicts with role, run, label, line, subset, out. comparisons: compare() outputs (selected vs
    run). c1: controls.compare_c1() outputs."""
    rows = []
    for r in records:
        for ms, vals in metric_sets(r["out"]).items():
            rows.append({"role": r["role"], "run": r["run"], "label": r.get("label", ""), "line": r["line"],
                         "split": r["out"]["split"], "subset": r.get("subset") or "", "metric_set": ms, **vals})
    for c in comparisons:
        jk = c["paired_jackknife"]
        rows.append({"role": "comparison", "run": c["incumbent"], "label": f"{c['candidate']} vs {c['incumbent']}",
                     "line": c["line"], "split": c["split"], "subset": "", "metric_set": "paired",
                     "dlogOR": c["delta_log_or"], "jk_CI95_lo": jk.get("ci95", [None, None])[0],
                     "jk_CI95_hi": jk.get("ci95", [None, None])[1], "boot_lb90": c["boot_delta_log_or_lb90"],
                     "dlogCMH": c["delta_log_cmh"], "boot_cmh_ub90": c["boot_delta_log_cmh_ub90"]})
    for c in c1:
        rows.append({"role": "control", "run": c["run"], "label": "C1 pair-specific gain", "line": c["line"],
                     "split": c["split"], "subset": "", "metric_set": "C1", "OR_run": c["or_run"],
                     "OR_C1": c["or_c1"], "dOR": c["delta_or"]["estimate"],
                     "dlogOR": c["delta_log_or"]["estimate"], "jk_CI95_lo": c["delta_log_or"]["ci95"][0],
                     "jk_CI95_hi": c["delta_log_or"]["ci95"][1], "c1_mode": c["c1_mode"]})
    return pd.DataFrame(rows)


def _fmt(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    if isinstance(v, (float, np.floating)):
        return f"{v:.3g}" if abs(v) < 1e4 else f"{v:.0f}"
    return str(v)


def to_markdown(df: pd.DataFrame, title: str = "Final evaluation") -> str:
    """One compact table per metric set (empty columns dropped)."""
    out = [f"# {title}", ""]
    key = ["role", "run", "label", "line", "split", "subset"]
    order = ["pooled", "cmh", "sign", "callset", "C1", "paired", "runtime"]
    sets = sorted(df["metric_set"].unique(), key=lambda s: (order.index(s) if s in order else len(order), s))
    for ms in sets:
        d = df[df["metric_set"] == ms].dropna(axis=1, how="all").drop(columns="metric_set")
        cols = [c for c in key if c in d.columns] + [c for c in d.columns if c not in key]
        out += [f"## {ms}", "", "| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
        out += ["| " + " | ".join(_fmt(v) for v in row) + " |" for row in d[cols].itertuples(index=False)]
        out.append("")
    return "\n".join(out)
