#!/usr/bin/env python
"""Profile the PIE X-Atlas/Orion benchmark from the compact caches (cache_pie_de.py) and the gene map
(build_gene_annotation.py). Writes data/processed/benchmark_profile.json and
data/processed/pie/target_summary_<cell>.tsv (per-target on-target knockdown and hit counts).

Definitions (explicit):
  tested pair   := (target, feature) row present in the DE table.
  on-target     := the feature whose identifier equals the target symbol, or, failing that, the feature
                   mapped to the same HGNC id (pie_gene_map.tsv).
  trans         := target and feature both placed on primary-assembly chromosomes and chrom differs.
                   (PAR genes use their Ensembl X location; unplaced ids are excluded from trans/cis.)
  significant   := p_adj < 0.05 or < 0.01 (BH as supplied by PIE; we verify below whether BH was
                   applied within target).
  DepMap-measurable feature := feature mapped (one-to-one) to a DepMap expression column.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pyarrow.dataset as ds
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data" / "processed"
PIE = PROC / "pie"
SRC = ROOT / "data" / "external" / "pie_x_atlas_orion" / "de"
CELLS = ["hct116", "hek293t"]
THR = {"0.05": -np.log10(0.05), "0.01": -np.log10(0.01)}


def q(x, qs=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)) -> dict:
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    return {str(k): float(v) for k, v in zip(qs, np.quantile(x, qs))} if len(x) else {}


def main() -> int:
    R: dict = {"definitions": __doc__.split("Definitions (explicit):")[1].strip()}
    tv = pd.read_csv(PIE / "vocab_targets.tsv", sep="\t", keep_default_na=False)
    fv = pd.read_csv(PIE / "vocab_features.tsv", sep="\t", keep_default_na=False)
    gm = pd.read_csv(PIE / ".." / "pie_gene_map.tsv", sep="\t", keep_default_na=False, na_values=[""],
                     dtype={"chrom": str}).set_index("pie_id")
    var = pd.read_csv(PROC / "depmap_cn_variability.tsv", sep="\t").set_index("depmap_column")

    t_chrom = tv.target.map(gm.chrom).to_numpy(dtype=object)
    f_chrom = fv.feature.map(gm.chrom).to_numpy(dtype=object)
    t_class = tv.target.map(gm.chrom_class)
    f_class = fv.feature.map(gm.chrom_class)
    t_placed = (t_class != "unplaced").to_numpy() & pd.notna(t_chrom)
    f_placed = (f_class != "unplaced").to_numpy() & pd.notna(f_chrom)
    chroms = sorted(set(x for x in t_chrom if isinstance(x, str)) | set(x for x in f_chrom if isinstance(x, str)))
    cid = {c: i for i, c in enumerate(chroms)}
    t_ci = np.array([cid.get(c, -1) if isinstance(c, str) else -1 for c in t_chrom], dtype=np.int16)
    f_ci = np.array([cid.get(c, -1) if isinstance(c, str) else -1 for c in f_chrom], dtype=np.int16)
    t_ci[~t_placed] = -1
    f_ci[~f_placed] = -1
    f_dep_expr = fv.feature.map(gm.depmap_expr_column).notna().to_numpy()
    t_dep_cn_col = tv.target.map(gm.depmap_cn_column)
    t_dep_cn = t_dep_cn_col.notna().to_numpy()
    # variable-CN flags for targets (two candidate definitions; see DATA_PROFILE)
    def var_flag(col, fn):
        if not isinstance(col, str):
            return False
        c = col.split(";")[0]
        return bool(c in var.index and fn(var.loc[c]))
    t_var_sd05 = np.array([var_flag(c, lambda r: r.sd_log2 > 0.5) for c in t_dep_cn_col])
    t_var_n1_20 = np.array([var_flag(c, lambda r: r["n_abs_gt_1"] >= 20) for c in t_dep_cn_col])
    t_var_any = np.array([var_flag(c, lambda r: r.frac_missing < 1) for c in t_dep_cn_col])

    # on-target feature code per target
    fcode = dict(zip(fv.feature, fv.code))
    hgnc_to_f = {}
    for f, code in zip(fv.feature, fv.code):
        h = gm.hgnc_id.get(f)
        if isinstance(h, str) and gm.map_method.get(f) != "collision":
            hgnc_to_f.setdefault(h, code)
    on_f = np.full(len(tv), -1, dtype=np.int32)
    on_how = np.array(["none"] * len(tv), dtype=object)
    for i, t in enumerate(tv.target):
        if t in fcode:
            on_f[i], on_how[i] = fcode[t], "identical_id"
        else:
            h = gm.hgnc_id.get(t)
            if isinstance(h, str) and h in hgnc_to_f:
                on_f[i], on_how[i] = hgnc_to_f[h], "same_hgnc"

    per_cell = {}
    target_sets = {}
    for cell in CELLS:
        T = pl.read_parquet(PIE / f"tested_{cell}.parquet")
        t = T["t"].to_numpy().astype(np.int64)
        f = T["f"].to_numpy().astype(np.int64)
        lfc = T["lfc"].to_numpy()
        mpadj = T["mlog10_padj"].to_numpy()
        n_t, n_f = len(tv), len(fv)
        tcount = np.bincount(t, minlength=n_t)
        targets = np.flatnonzero(tcount > 0)
        target_sets[cell] = set(tv.target.iloc[targets])
        fcount = np.bincount(f, minlength=n_f)
        trans = (t_ci[t] >= 0) & (f_ci[f] >= 0) & (t_ci[t] != f_ci[f])
        cis = (t_ci[t] >= 0) & (t_ci[t] == f_ci[f])
        is_on = on_f[t] == f
        res = {"n_rows": int(len(t)), "n_targets": int(len(targets)), "n_features": int((fcount > 0).sum()),
               "features_per_target": q(tcount[targets]),
               "targets_per_feature": q(fcount[fcount > 0]),
               "frac_pairs_trans": float(trans.mean()), "frac_pairs_cis": float(cis.mean()),
               "nan_lfc": int(np.isnan(lfc).sum()), "nan_padj": int(np.isnan(mpadj).sum())}
        tinfo = pd.read_csv(PIE / f"targets_{cell}.tsv", sep="\t", keep_default_na=False).set_index("target")
        res["target_ncells"] = q(tinfo.target_ncells)
        res["ref_ncells"] = int(tinfo.ref_ncells.max())
        finfo = pd.read_csv(PIE / f"features_{cell}.tsv", sep="\t", keep_default_na=False).set_index("feature")
        res["spearman_feature_ref_mean_vs_n_targets_tested"] = float(
            finfo.ref_mean.rank().corr(finfo.n_targets.rank()))
        res["frac_features_tested_for_all_targets"] = float((finfo.n_targets == len(targets)).mean())

        # on-target knockdown
        sel = is_on
        kd = pd.DataFrame({"t": t[sel], "lfc": lfc[sel], "mpadj": mpadj[sel]}).set_index("t")
        on_lfc = kd.lfc.reindex(targets).to_numpy()
        on_p = kd.mpadj.reindex(targets).to_numpy()
        has = ~np.isnan(on_lfc)
        res["on_target"] = {
            "n_targets_with_own_feature_tested": int(has.sum()),
            "n_own_feature_absent": int((~has).sum()),
            "match_how": {k: int(v) for k, v in pd.Series(on_how[targets][has]).value_counts().items()},
            "lfc_q": q(on_lfc[has]),
            "frac_lfc_lt_0": float((on_lfc[has] < 0).mean()),
            "frac_lfc_lt_-1": float((on_lfc[has] < -1).mean()),
            "frac_padj_lt_0.05_and_neg": float(((on_p[has] > THR["0.05"]) & (on_lfc[has] < 0)).mean()),
            "frac_padj_lt_0.05_and_pos": float(((on_p[has] > THR["0.05"]) & (on_lfc[has] > 0)).mean()),
            "frac_padj_lt_0.01_and_neg": float(((on_p[has] > THR["0.01"]) & (on_lfc[has] < 0)).mean()),
            "frac_lfc_lt_-20_(perturbed_mean_zero;pseudocount_artefact)": float((on_lfc[has] < -20).mean()),
        }
        res["lfc_extremes"] = {"frac_abs_lfc_gt_10": float((np.abs(lfc) > 10).mean()),
                               "frac_abs_lfc_gt_10_among_padj<0.05": float((np.abs(lfc[mpadj > THR["0.05"]]) > 10).mean())}
        # "consistently tested" features: tested for every target of the cell line
        f_all = fcount == len(targets)
        res["n_features_tested_all_targets"] = int(f_all.sum())
        res["n_features_tested_all_targets_in_depmap_expr"] = int((f_all & f_dep_expr).sum())
        res["features_partial_ref_mean_max"] = float(finfo.loc[[x for x in fv.feature[(fcount > 0) & ~f_all]], "ref_mean"].max())
        res["features_all_ref_mean_min"] = float(finfo.loc[[x for x in fv.feature[f_all]], "ref_mean"].min())

        # hits per target
        summ = pd.DataFrame({"target": tv.target.iloc[targets].to_numpy(), "n_tested": tcount[targets],
                             "on_target_lfc": on_lfc, "on_target_mlog10_padj": on_p,
                             "target_ncells": tinfo.target_ncells.reindex(tv.target.iloc[targets]).to_numpy(),
                             "chrom": t_chrom[targets], "in_depmap_cn": t_dep_cn[targets],
                             "cn_var_sd_gt_0.5": t_var_sd05[targets], "cn_var_n_abs1_ge20": t_var_n1_20[targets]})
        hits = {}
        for lab, thr in THR.items():
            sig = mpadj > thr
            for scope, msk in [("all", np.ones_like(sig)), ("trans", trans), ("cis", cis),
                               ("trans_depmap_expr", trans & f_dep_expr[f])]:
                s_ = sig & msk & ~is_on
                cnt = np.bincount(t[s_], minlength=n_t)[targets]
                up = np.bincount(t[s_ & (lfc > 0)], minlength=n_t)[targets]
                summ[f"n_sig{lab}_{scope}"] = cnt
                if scope in ("all", "trans"):
                    summ[f"n_up{lab}_{scope}"] = up
                hits[f"padj<{lab}|{scope}"] = {
                    "total_pairs": int(cnt.sum()), "targets_with_ge1": int((cnt >= 1).sum()),
                    "targets_with_ge10": int((cnt >= 10).sum()), "targets_with_ge100": int((cnt >= 100).sum()),
                    "per_target_q": q(cnt, (0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)),
                    "mean_per_target": float(cnt.mean()),
                    "frac_up": float(up.sum() / max(cnt.sum(), 1)),
                }
        res["hits_excluding_on_target"] = hits
        # targets whose own gene is expressed (own feature tested) vs not (CRISPRi of a non-expressed gene
        # = natural negative control)
        sig_tr = (mpadj > THR["0.05"]) & trans & ~is_on
        cnt = np.bincount(t[sig_tr], minlength=n_t)[targets]
        cnt_c = np.bincount(t[sig_tr & f_all[f]], minlength=n_t)[targets]
        summ["n_sig0.05_trans_consistent_features"] = cnt_c
        expd = has
        res["trans_hits_by_target_expression"] = {
            k: {"n_targets": int(msk.sum()), "mean": float(cnt[msk].mean()), "median": float(np.median(cnt[msk])),
                "frac_ge1": float((cnt[msk] >= 1).mean()), "frac_ge10": float((cnt[msk] >= 10).mean()),
                "mean_consistent_features": float(cnt_c[msk].mean()), "total": int(cnt[msk].sum())}
            for k, msk in [("own_gene_expressed", expd), ("own_gene_not_tested", ~expd)]}
        for k, msk in [("own_gene_expressed", expd), ("own_gene_not_tested", ~expd)]:
            up = np.bincount(t[sig_tr & (lfc > 0)], minlength=n_t)[targets]
            res["trans_hits_by_target_expression"][k]["frac_up"] = float(up[msk].sum() / max(cnt[msk].sum(), 1))
        res["frac_trans_hits_in_partially_tested_features"] = float(
            (sig_tr & ~f_all[f]).sum() / max(sig_tr.sum(), 1))
        # concentration: share of trans hits from top 1% / 10% of targets
        c = np.sort(summ["n_sig0.05_trans"].to_numpy())[::-1]
        res["trans_hit_concentration"] = {"top1pct_share": float(c[: max(1, len(c) // 100)].sum() / max(c.sum(), 1)),
                                          "top10pct_share": float(c[: max(1, len(c) // 10)].sum() / max(c.sum(), 1))}
        # knockdown strength vs hit count
        ok = summ.on_target_lfc.notna()
        res["spearman_on_target_lfc_vs_n_trans_hits"] = float(
            summ.loc[ok, "on_target_lfc"].rank().corr(summ.loc[ok, "n_sig0.05_trans"].rank()))
        res["spearman_target_ncells_vs_n_trans_hits"] = float(
            summ.target_ncells.rank().corr(summ["n_sig0.05_trans"].rank()))
        summ.to_csv(PIE / f"target_summary_{cell}.tsv", sep="\t", index=False)

        # universe preview (NOT the final benchmark definition): targets with DepMap CN x trans features
        # measured in DepMap expression AND tested for every target in this line, restricted by candidate
        # variable-CN definitions
        t_expr = np.zeros(n_t, dtype=bool)
        t_expr[targets[has]] = True
        uni = {}
        for name, tmask in [("target_in_depmap_cn", t_dep_cn & t_var_any), ("cn_sd_gt_0.5", t_var_sd05),
                            ("cn_n_abs1_ge20", t_var_n1_20),
                            ("target_in_depmap_cn&own_gene_expressed", t_dep_cn & t_var_any & t_expr),
                            ("cn_n_abs1_ge20&own_gene_expressed", t_var_n1_20 & t_expr)]:
            msk = tmask[t] & trans & f_dep_expr[f] & ~is_on & f_all[f]
            pos05 = msk & (mpadj > THR["0.05"])
            same_sign = pos05 & (lfc < 0)
            uni[name] = {"n_targets": int(len(np.unique(t[msk]))), "n_pairs": int(msk.sum()),
                         "n_pos_padj<0.05": int(pos05.sum()), "pos_rate": float(pos05.sum() / max(msk.sum(), 1)),
                         "frac_pos_down": float(same_sign.sum() / max(pos05.sum(), 1))}
        res["universe_preview"] = uni
        per_cell[cell] = res
        del T, t, f, lfc, mpadj

    # target overlap
    a, b = target_sets["hct116"], target_sets["hek293t"]
    R["target_overlap"] = {"hct116": len(a), "hek293t": len(b), "both": len(a & b), "hct116_only": len(a - b),
                           "hek293t_only": len(b - a), "union": len(a | b)}
    both = sorted(a & b)
    # cross-line reproducibility of on-target KD and hit counts
    s1 = pd.read_csv(PIE / "target_summary_hct116.tsv", sep="\t", keep_default_na=False).set_index("target")
    s2 = pd.read_csv(PIE / "target_summary_hek293t.tsv", sep="\t", keep_default_na=False).set_index("target")
    for c in ["on_target_lfc", "n_sig0.05_trans"]:
        x1 = pd.to_numeric(s1.loc[both, c], errors="coerce")
        x2 = pd.to_numeric(s2.loc[both, c], errors="coerce")
        R["target_overlap"][f"spearman_{c}_across_lines"] = float(x1.rank().corr(x2.rank()))
    R["per_cell"] = per_cell

    # BH scope check: recompute BH within target for a few targets from the raw parquet
    bh = {}
    for cell in CELLS:
        summ = pd.read_csv(PIE / f"target_summary_{cell}.tsv", sep="\t", keep_default_na=False)
        picks = summ.sort_values("n_sig0.05_all").target.iloc[[len(summ) // 2, -1]].tolist()
        dset = ds.dataset(SRC / f"{cell}.parquet")
        out = {}
        for tg in picks:
            tab = dset.to_table(columns=["p_value", "p_adj"], filter=ds.field("target") == tg).to_pandas()
            within = multipletests(tab.p_value, method="fdr_bh")[1]
            out[tg] = {"n": len(tab), "max_abs_diff_within_target_BH": float(np.max(np.abs(within - tab.p_adj)))}
        bh[cell] = out
    R["bh_scope_check"] = bh

    # identifier / DepMap overlap
    ov = {}
    for cell in CELLS:
        tt = gm[gm[f"target_{cell}"]]
        ff = gm[gm[f"feature_{cell}"]]
        ov[cell] = {"targets": len(tt), "targets_in_depmap_cn": int(tt.depmap_cn_column.notna().sum()),
                    "targets_in_depmap_expr": int(tt.depmap_expr_column.notna().sum()),
                    "targets_in_depmap_cn_and_expr": int((tt.depmap_cn_column.notna() & tt.depmap_expr_column.notna()).sum()),
                    "targets_cn_sd_gt_0.5": int(t_var_sd05[tv.target.isin(tt.index).to_numpy()].sum()),
                    "targets_cn_n_abs1_ge20": int(t_var_n1_20[tv.target.isin(tt.index).to_numpy()].sum()),
                    "features": len(ff), "features_in_depmap_expr": int(ff.depmap_expr_column.notna().sum()),
                    "features_protein_coding": int((ff.locus_group == "protein-coding gene").sum())}
    R["depmap_overlap"] = ov

    with open(PROC / "benchmark_profile.json", "w") as fh:
        json.dump(R, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in R.items() if k != "definitions"}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
