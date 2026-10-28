#!/usr/bin/env python
"""Worked example: trace ONE HCT116 TRAIN regulator through the selected CN-only method.

    .venv/bin/python docs/explainer/make_worked_example.py            # -> docs/explainer/worked_example.json
    .venv/bin/python docs/explainer/make_worked_example.py --regulator "PPP1R15B (84919)"

Selected method (results/final/selected.yaml): g1_t100_w0__a3_r095_s2_w0 = M1 (LOCO reduced form) + G1 (t_cis >= 10)
+ A3 (SuSiE L=1, rho 0.95, s 2, w_co 0). Every intermediate is recomputed here with the production functions in
src/depmaptr/ and checked against the committed score run data/results/scores/g1_t100_w0__a3_r095_s2_w0/scores.npz.

Reads only: data/processed/{panel_s2,frozen_flags,gene_annotation}.tsv, the matched CN/expression parquets,
data/processed/bench/{universe_*_hct116.tsv, positives_hct116_TRAIN.parquet} (TRAIN labels only), and the score
files of the selected run and of the context layer (…__i1_fuse_g2). No TEST file is opened; the script refuses a
regulator that is not a HCT116 TRAIN regulator. Nothing in the repository is modified except the output JSON.
"""
from __future__ import annotations

import argparse
import importlib.util
import warnings
import json
import struct
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

# macOS Accelerate BLAS emits spurious "divide by zero / overflow in matmul" RuntimeWarnings on finite inputs;
# every result below is checked against statsmodels and the committed score file instead.
warnings.filterwarnings("ignore", message=".*encountered in matmul", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from depmaptr import attribution, gates  # noqa: E402
from depmaptr import io as dio  # noqa: E402
from depmaptr.assoc import Residualiser, centred_norm, diag_partial  # noqa: E402
from depmaptr.covariates import DesignBuilder  # noqa: E402
from depmaptr.pipeline import GroupCache, Params, group_t, prepare, target_context_percentile  # noqa: E402

SEL = "g1_t100_w0__a3_r095_s2_w0"
CTX = "g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2"
SCORES = ROOT / "data" / "results" / "scores"
BENCH = ROOT / "data" / "processed" / "bench"
HCT116 = "ACH-000971"


# ----------------------------------------------------------------------------- helpers
def mmap_npz(path: Path, key: str) -> np.ndarray:
    """Memory-map one member of an uncompressed .npz (np.savez) without loading it."""
    z = zipfile.ZipFile(path)
    info = z.getinfo(key + ".npy")
    assert info.compress_type == zipfile.ZIP_STORED, f"{key} is compressed"
    with open(path, "rb") as fh:
        fh.seek(info.header_offset)
        h = fh.read(30)
        n, m = struct.unpack("<HH", h[26:30])
        fh.seek(info.header_offset + 30 + n + m)
        ver = np.lib.format.read_magic(fh)
        shape, fortran, dtype = np.lib.format._read_array_header(fh, ver)
        off = fh.tell()
    return np.memmap(path, dtype=dtype, mode="r", offset=off, shape=shape, order="F" if fortran else "C")


def r4(x, nd=4):
    """Round floats (recursively) to nd significant digits for a compact JSON."""
    if isinstance(x, dict):
        return {k: r4(v, nd) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [r4(v, nd) for v in x]
    if isinstance(x, (np.floating, float)):
        x = float(x)
        if not np.isfinite(x):
            return None
        return float(f"{x:.{nd}g}")
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def sym(key: str) -> str:
    return key.split(" (")[0]


def selected_params() -> Params:
    spec = importlib.util.spec_from_file_location("run_model", ROOT / "scripts" / "run_model.py")
    rm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rm)
    defaults, table = rm.load_variants(ROOT / "results" / "variants.yaml")
    kind, params, _ = rm.resolve("g1_t100_w0+a3_r095_s2_w0", defaults, table)
    assert kind == "model"
    return Params.from_dict(params)


def ols_t(y, cols):
    X = np.column_stack(cols)
    fit = sm.OLS(y, X).fit()
    return float(fit.tvalues[-1]), float(fit.params[-1]), int(fit.df_resid)


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--regulator", default="PPP1R15B (84919)")
    ap.add_argument("--targets", default="MYBL2 (4605),NT5E (4907)",
                    help="named targets for the regression walk-through (a confounded one is added automatically)")
    ap.add_argument("--out", default=str(Path(__file__).with_name("worked_example.json")))
    a = ap.parse_args(argv)
    REG = a.regulator
    timing = {}
    T0 = time.time()

    # ---- benchmark context (TRAIN only) -------------------------------------------------------------------
    ureg = pd.read_csv(BENCH / "universe_regulators_hct116.tsv", sep="\t", dtype={"chrom": str})
    utgt = pd.read_csv(BENCH / "universe_targets_hct116.tsv", sep="\t", dtype={"chrom": str})
    row = ureg.loc[ureg.regulator == REG]
    if row.empty or row.split.iloc[0] != "TRAIN":
        raise SystemExit(f"{REG} is not a HCT116 TRAIN universe regulator; refusing (TEST labels are sealed)")
    train_regs = ureg.loc[ureg.split == "TRAIN"].reset_index(drop=True)
    pos = pd.read_parquet(BENCH / "positives_hct116_TRAIN.parquet")
    pos = pos[pos.pos_primary]
    pos_reg = pos[pos.regulator == REG].set_index("target")

    # ---- data + design --------------------------------------------------------------------------------------
    t = time.time()
    ds = dio.load_real("S2", checksums=False)
    timing["load_real_s"] = time.time() - t
    p = selected_params()
    t = time.time()
    prep = prepare(ds, p)
    timing["prepare_s"] = time.time() - t
    genes = ds.genes
    c = str(genes.loc[REG, "chrom"])
    panel = ds.panel
    names = panel["CellLineName"].astype(str).to_numpy()
    lineage = panel["lineage"].astype(str).to_numpy()
    i_hct = int(np.flatnonzero(panel["ModelID"].to_numpy() == HCT116)[0])

    # ---- step 1: the regulator's CN across the panel -----------------------------------------------------
    jA = prep.cn_idx[REG]
    relcn = ds.relcn[:, jA]
    xA = prep.x_all[:, jA]
    yA = prep.y_all[:, prep.ex_idx[REG]]
    qs = [0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1]
    order = np.argsort(xA)
    ex = lambda idx: [{"line": names[k], "lineage": lineage[k], "relCN": float(relcn[k]), "x": float(xA[k]),
                       "expr": float(yA[k])} for k in idx]
    hist_edges = np.arange(-2.0, 3.01, 0.25)
    cn_block = {
        "n_lines": int(ds.n),
        "transform": "x = log2(clip(relCN, 2^-2, 2^3)); relCN = 1 is copy-neutral relative to the line's ploidy",
        "quantiles_x": {f"q{int(q * 100):02d}": float(np.quantile(xA, q)) for q in qs},
        "quantiles_relCN": {f"q{int(q * 100):02d}": float(np.quantile(relcn, q)) for q in qs},
        "n_gain_x_ge_1": int((xA >= 1).sum()), "n_loss_x_le_-1": int((xA <= -1).sum()),
        "n_abs_x_ge_1 (variable_cn statistic; >=20 required)": int((np.abs(xA) >= 1).sum()),
        "n_at_floor_relCN_lt_0.25": int((relcn < 0.25).sum()),
        "hct116": ex([i_hct])[0],
        "lowest_3": ex(order[:3]), "highest_3": ex(order[-3:][::-1]),
        "hist": {"edges": hist_edges.tolist(), "counts": np.histogram(xA, bins=hist_edges)[0].tolist()},
    }

    # ---- step 2: cis instrument (own expression vs own CN) -> G1 -----------------------------------------
    t = time.time()
    gv = GroupCache(prep).get(c)        # D_c residualiser + residualised targets/cis basis for chromosome c
    timing["residualise_group_s"] = time.time() - t
    res = gv["res"]
    D_c, Dnames = prep.builder.design(c)
    xr, yr = res(xA[:, None]), res(yA[:, None])
    r_cis, t_cis, df_cis = diag_partial(xr, yr, df_base=res.df_base, u_raw=centred_norm(xA[:, None]),
                                        y_raw=centred_norm(yA[:, None]))
    t_cis = float(t_cis[0])
    raw_r = float(np.corrcoef(xA, yA)[0, 1])
    raw_slope = float(np.polyfit(xA, yA, 1)[0])
    adj_slope = float((xr[:, 0] @ yr[:, 0]) / (xr[:, 0] @ xr[:, 0]))
    sm_t, sm_b, sm_df = ols_t(yA, [D_c, xA])
    g1 = bool(gates.g1_pass(np.array([t_cis]), p.g1_tau)[0])
    s_A = int(gates.s_cis(np.array([t_cis]))[0])
    keep = np.r_[order[:: max(1, ds.n // 400)], i_hct]          # ~400-point subsample for plots
    cis_block = {
        "model": "y_A ~ x_A | D_c (D_c = intercept + 56 shared covariates + 2 LOCO burden + arm CN of every arm not on A's chromosome)",
        "design_columns": int(D_c.shape[1]), "design_rank": int(res.rank), "df": float(df_cis),
        "raw_pearson_r": raw_r, "raw_slope_log2TPM_per_log2CN": raw_slope,
        "adjusted_partial_r": float(r_cis[0]), "adjusted_slope": adj_slope, "t_cis": t_cis, "F_first_stage": t_cis ** 2,
        "statsmodels_check": {"t": sm_t, "slope": sm_b, "df_resid": sm_df},
        "G1_rule": "pass iff t_cis > 0 and t_cis^2 >= 100 (tau = 100, i.e. t_cis >= 10); fail -> weight 0",
        "G1_pass": g1, "s_A": s_A, "s_A_rule": "s_A = -1 iff t_cis <= -3, else +1",
        "scatter_subsample": {"x": [round(float(v), 3) for v in xA[keep]], "y": [round(float(v), 3) for v in yA[keep]],
                              "hct116_last": True},
    }

    # ---- step 3: A3 segmentation of A's chromosome -------------------------------------------------------
    t = time.time()
    chrom_keys = np.asarray([k for k in ds.cn_keys if str(genes.loc[k, "chrom"]) == c])
    loco = DesignBuilder(ds, prep.Z, prep.arms, "loco")
    Dl, _ = loco.design(c)
    xres = Residualiser(Dl)(prep.x_all[:, [prep.cn_idx[k] for k in chrom_keys]])
    labels = attribution.greedy_segments(xres, p.a3_rho)
    timing["segment_s"] = time.time() - t
    pos_in = int(np.flatnonzero(chrom_keys == REG)[0])
    seg = labels[pos_in]
    members = chrom_keys[labels == seg]

    # ---- step 4: reduced-form t for every segment member x every target (production path) ----------------
    t = time.time()
    tS, dfS, _ = group_t(prep, gv, members, None)
    timing["trans_t_segment_s"] = time.time() - t
    cisS = np.zeros(len(members))
    has_e = np.array([k in prep.ex_idx for k in members])
    if has_e.any():
        ke = members[has_e]
        xe = prep.x_all[:, [prep.cn_idx[k] for k in ke]]
        ye = prep.y_all[:, [prep.ex_idx[k] for k in ke]]
        _, tc, _ = diag_partial(res(xe), res(ye), df_base=res.df_base, u_raw=centred_norm(xe), y_raw=centred_norm(ye))
        cisS[has_e] = np.nan_to_num(tc, nan=0.0)
    alphaS = attribution.segment_alpha(np.zeros(len(members), int), tS, dfS, cisS, p.a3_s, sigma0=p.a3_sigma0)
    iA = int(np.flatnonzero(members == REG)[0])

    targets = np.asarray(ds.targets)
    tch = genes.loc[targets, "chrom"].astype(str).to_numpy()
    trans = tch != c
    tA, aA = tS[iA].copy(), alphaS[iA]
    tA[targets == REG] = np.nan
    w1 = gates.g1_weight(np.array([t_cis]), p.g1_tau, p.g1_w)[0]
    S = (tA * w1 * aA).astype(np.float32)

    # ---- check against the committed score file ----------------------------------------------------------
    f_sel = SCORES / SEL / "scores.npz"
    z = np.load(f_sel)
    regs_all = z["regulators"]
    assert (z["targets"] == targets).all()
    ri = int(np.flatnonzero(regs_all == REG)[0])
    Ssel, Tsel, Asel = mmap_npz(f_sel, "score"), mmap_npz(f_sel, "t"), mmap_npz(f_sel, "alpha")
    s_file, t_file, a_file = np.asarray(Ssel[ri]), np.asarray(Tsel[ri]), np.asarray(Asel[ri])
    fin = np.isfinite(s_file) & trans
    check = {
        "max_abs_diff_t": float(np.nanmax(np.abs(t_file[fin] - tA[fin]))),
        "max_abs_diff_alpha": float(np.max(np.abs(a_file[fin] - aA[fin]))),
        "max_abs_diff_score": float(np.max(np.abs(s_file[fin] - S[fin]))),
        "cis_t_file": float(z["cis_t"][ri]), "s_cis_file": int(z["s_cis"][ri]),
    }
    reg_set = set(ds.regulators)
    univ_regs = set(ureg.regulator)
    split_of = ureg.set_index("regulator")["split"].to_dict()
    gmap = pd.read_csv(ROOT / "data" / "processed" / "pie_gene_map.tsv", sep="\t", low_memory=False)
    perturbed = set(gmap.loc[gmap["target_hct116"].astype(str) == "True", "depmap_cn_column"].dropna())
    seg_rows = []
    for k, g in enumerate(members):
        g1g = bool(gates.g1_pass(np.array([cisS[k]]), p.g1_tau)[0])
        seg_rows.append({
            "gene": g, "start_bp": int(genes.loc[g, "start"]) if pd.notna(genes.loc[g, "start"]) else None,
            "has_expression": bool(has_e[k]), "t_cis": float(cisS[k]), "G1_pass": g1g,
            "log_prior_unnorm (t_cis/s)": float(cisS[k] / p.a3_s),
            "prior_pi": float(np.exp(cisS[k] / p.a3_s - np.logaddexp.reduce(cisS / p.a3_s))),
            "alpha_median_trans": float(np.median(alphaS[k][trans])),
            "alpha_q10_q90": [float(np.quantile(alphaS[k][trans], 0.1)), float(np.quantile(alphaS[k][trans], 0.9))],
            "is_tool_regulator (variable_cn)": g in reg_set,
            "perturbed_in_HCT116_screen": g in perturbed,
            "hct116_benchmark": ("TRAIN" if split_of.get(g) == "TRAIN" else ("held-out" if g in univ_regs else "not in universe")),
            "TRAIN_hits (p_adj<0.05)": int((pos.regulator == g).sum()) if split_of.get(g) == "TRAIN" else None,
            "r_CN_with_regulator_resid": float(np.corrcoef(xres[:, pos_in], xres[:, np.flatnonzero(chrom_keys == g)[0]])[0, 1]),
        })
    # contrast neighbour: the member (other than A) with the lowest median alpha
    others = [k for k in range(len(members)) if k != iA]
    k_low = min(others, key=lambda k: seg_rows[k]["alpha_median_trans"]) if others else None

    # ---- step 5: walk-through targets --------------------------------------------------------------------
    jT = {k: j for j, k in enumerate(targets)}
    univ_t = set(utgt.target)
    named = [s for s in a.targets.split(",") if s]
    # a confounded example: the trans universe target with the largest |raw r| whose adjusted |t| < 2
    YT = prep.YT
    raw_r_all = np.full(len(targets), np.nan)
    xc = (xA - xA.mean()) / xA.std()
    Yc = (YT - YT.mean(0)) / YT.std(0)
    raw_r_all = (xc @ Yc) / ds.n
    cand = np.flatnonzero(trans & np.isfinite(tA) & (np.abs(tA) < 2) & np.isin(targets, list(univ_t)))
    jconf = cand[np.argmax(np.abs(raw_r_all[cand]))]
    named.append(targets[jconf])
    arm_cols = [i for i, nm in enumerate(Dnames) if nm.startswith("arm[")]
    D_noarm = np.delete(D_c, arm_cols, axis=1)
    walk = []
    for B in named:
        j = jT[B]
        yB, xB, dB = prep.YT[:, j], prep.XT[:, j], prep.DEL[:, j]
        W = [xB] + ([dB] if dB.any() else [])
        t_full, b_full, df_full = ols_t(yB, [D_c] + W + [xA])
        t_noarm, _, _ = ols_t(yB, [D_noarm] + W + [xA])
        t_nocov, _, _ = ols_t(yB, [np.ones(ds.n)] + W + [xA])
        rr = float(raw_r_all[j])
        hit = pos_reg.lfc_sign.get(B) if B in pos_reg.index else None
        kd = int(-np.sign(S[j]) * s_A)
        walk.append({
            "target": B, "chrom": tch[j], "in_hct116_universe": B in univ_t,
            "raw_pearson_r": rr, "raw_t (naive, df=n-2; baseline B1)": rr * np.sqrt((ds.n - 2) / (1 - rr * rr)),
            "t_cisB_only (y_B ~ 1 + x_B [+del_B] + x_A)": t_nocov,
            "t_no_arm_block (all covariates except arm CN)": t_noarm,
            "t_full_statsmodels": t_full, "beta_full_log2TPM_per_log2CN": b_full, "df_resid": df_full,
            "t_engine": float(tA[j]), "abs_diff_engine_vs_statsmodels": abs(float(tA[j]) - t_full),
            "homdel_indicator_in_model": bool(dB.any()),
            "alpha": float(aA[j]), "S": float(S[j]), "pred_kd_sign": kd,
            "TRAIN_hit": hit is not None, "observed_lfc_sign": None if hit is None else int(hit),
            "concordant": None if hit is None else bool(kd == int(hit)),
            "alpha_by_member": [{"gene": members[k], "t": float(tS[k, j]), "df": float(dfS[k, j]),
                                 "log_ABF": float(attribution.log_abf(np.nan_to_num(tS[k, j]), dfS[k, j], p.a3_sigma0)),
                                 "log_prior": float(cisS[k] / p.a3_s), "alpha": float(alphaS[k, j])}
                                for k in range(len(members))],
        })

    # residual scatter for the first named target (the adjusted association, after FWL)
    j0 = jT[named[0]]
    yb0 = gv["Y"][:, j0]
    Wb = np.column_stack([gv["W"][0][:, j0]] + ([gv["W"][1][:, j0]] if prep.DEL[:, j0].any() else []))
    Qb, _ = np.linalg.qr(Wb)
    x_part = xr[:, 0] - Qb @ (Qb.T @ xr[:, 0])
    y_part = yb0 - Qb @ (Qb.T @ yb0)
    ks = order[:: max(1, ds.n // 400)]
    fwl_scatter = {"target": named[0], "raw_x": [round(float(v), 3) for v in xA[ks]],
                   "raw_y": [round(float(v), 3) for v in prep.YT[ks, j0]],
                   "adj_x": [round(float(v), 3) for v in x_part[ks]], "adj_y": [round(float(v), 3) for v in y_part[ks]]}

    # ---- step 6: top targets and TRAIN labels ------------------------------------------------------------
    def top_table(mask, n=10):
        sc = np.where(mask & np.isfinite(S), np.abs(S), -1)
        out = []
        for j in np.argsort(-sc)[:n]:
            B = targets[j]
            hit = pos_reg.lfc_sign.get(B) if B in pos_reg.index else None
            kd = int(-np.sign(S[j]) * s_A)
            out.append({"target": B, "chrom": tch[j], "t": float(tA[j]), "alpha": float(aA[j]), "S": float(S[j]),
                        "pred_kd": "down" if kd < 0 else "up", "in_hct116_universe": B in univ_t,
                        "TRAIN_hit": hit is not None,
                        "observed": None if hit is None else ("down" if hit < 0 else "up"),
                        "concordant": None if hit is None else bool(kd == int(hit))})
        return out

    in_univ = np.isin(targets, list(univ_t))
    top_all = top_table(trans)
    top_univ = top_table(trans & in_univ)

    # TRAIN top-K (K = 0.1% of TRAIN universe pairs) from the committed scores, labels used only to count hits
    ri_tr = np.array([int(np.flatnonzero(regs_all == r)[0]) for r in train_regs.regulator])
    tj_u = np.array([jT[k] for k in utgt.target])
    M = np.abs(np.asarray(Ssel[ri_tr][:, tj_u], dtype=np.float32))
    same = train_regs.chrom.to_numpy()[:, None] == utgt.chrom.to_numpy()[None, :]
    M[same] = -np.inf
    M[np.isnan(M)] = -1.0          # NaN -> minimum (harness rule)
    Npairs = int((~same).sum())
    K = int(round(0.001 * Npairs))
    thr = float(np.partition(M.ravel(), -K)[-K])
    n_at_thr = int((M >= thr).sum())
    ntop = (M >= thr).sum(1)
    rk = int(np.flatnonzero(train_regs.regulator.to_numpy() == REG)[0])
    rank_by_ntop = int((ntop > ntop[rk]).sum() + 1)
    in_top = M[rk] >= thr
    tgt_u = utgt.target.to_numpy()
    trans_u = ~same[rk]
    is_hit = np.isin(tgt_u, pos_reg.index.to_numpy())
    a_ = int((in_top & is_hit).sum()); b_ = int((in_top & ~is_hit).sum())
    c_ = int((~in_top & is_hit & trans_u).sum()); d_ = int((~in_top & ~is_hit & trans_u).sum())
    conc = []
    for j in np.flatnonzero(in_top & is_hit):
        B = tgt_u[j]
        conc.append(int(-np.sign(S[jT[B]]) * s_A) == int(pos_reg.lfc_sign[B]))
    pos_all = pos.groupby("regulator").size()
    deg_rank = int((pos_all > pos_all.get(REG, 0)).sum() + 1)
    # label-free genome-wide rank: regulators ordered by max |S| over trans targets
    regchrom = genes.loc[regs_all, "chrom"].astype(str).to_numpy()
    maxabs = np.zeros(len(regs_all), np.float32)
    for cc in pd.unique(regchrom):
        rows = np.flatnonzero(regchrom == cc)
        blk = np.abs(np.asarray(Ssel[rows[0]:rows[-1] + 1]))[:, tch != cc]
        maxabs[rows] = np.nanmax(np.where(np.isfinite(blk), blk, 0), axis=1)
    rank_block = {
        "train_universe_pairs": Npairs, "K_top_0.1pct": K, "abs_S_threshold": thr, "pairs_at_or_above_threshold": n_at_thr,
        "regulator_pairs_in_topK": int(in_top.sum()), "rank_among_3806_TRAIN_regulators_by_topK_pairs": rank_by_ntop,
        "TRAIN_regulators_with_any_topK_pair": int((ntop > 0).sum()),
        "regulator_TRAIN_hits_total": int(len(pos_reg)), "regulator_rank_by_TRAIN_hit_count": deg_rank,
        "trans_universe_targets_for_regulator": int(trans_u.sum()),
        "within_regulator_2x2 [topK&hit, topK&nonhit, other&hit, other&nonhit]": [a_, b_, c_, d_],
        "hit_rate_in_topK": a_ / max(1, a_ + b_), "hit_rate_outside_topK": c_ / max(1, c_ + d_),
        "within_regulator_OR": (a_ * d_) / max(1, b_ * c_),
        "sign_concordant_among_topK_hits": int(sum(conc)), "n_topK_hits": len(conc),
        "TRAIN_wide_hit_rate": float(len(pos) / Npairs),
        "genome_wide_rank_by_max_abs_S (of 12141 regulators; label-free)": int((maxabs > maxabs[ri]).sum() + 1),
        "n_regulators_passing_G1_genome_wide": int((z["cis_t"] >= 10).sum()),
        "n_TRAIN_regulators_passing_G1": int((z["cis_t"][ri_tr] >= 10).sum()),
    }
    del M

    # ---- step 7: context layer (reported variant, not the selected method) -------------------------------
    f_ctx = SCORES / CTX / "scores.npz"
    s_ctx = np.asarray(mmap_npz(f_ctx, "score")[int(np.flatnonzero(np.load(f_ctx)["regulators"] == REG)[0])])
    rho_T = target_context_percentile(ds, "HCT116")
    rho_P_back = np.abs(s_ctx) / rho_T ** 2
    # exact average-rank percentile of |S| among all finite trans entries (incl. zeros) for 3 targets
    q_targets = [top_univ[0]["target"], named[0], named[1]]
    qv = np.array([abs(S[jT[B]]) for B in q_targets], np.float64)
    n_less = np.zeros(len(qv)); n_eq = np.zeros(len(qv)); N = 0; n_zero = 0
    for cc in pd.unique(regchrom):
        rows = np.flatnonzero(regchrom == cc)
        for s0 in range(rows[0], rows[-1] + 1, 1024):
            blk = np.abs(np.asarray(Ssel[s0:min(rows[-1] + 1, s0 + 1024)]))[:, tch != cc].ravel()
            blk = blk[np.isfinite(blk)]
            N += blk.size
            n_zero += int(np.count_nonzero(blk == 0))
            for qi, v in enumerate(qv.astype(np.float32)):
                n_less[qi] += np.count_nonzero(blk < v)
                n_eq[qi] += np.count_nonzero(blk == v)
    rho_P_exact = (n_less + (n_eq + 1) / 2) / N
    n_greater = N - n_less - n_eq
    ctx_rows = []
    sc = np.where(trans & in_univ & np.isfinite(s_ctx), np.abs(s_ctx), -1)
    for j in np.argsort(-sc)[:10]:
        B = targets[j]
        hit = pos_reg.lfc_sign.get(B) if B in pos_reg.index else None
        ctx_rows.append({"target": B, "S_ctx": float(s_ctx[j]), "rho_T(HCT116 expr pct)": float(rho_T[j]),
                         "rho_P": float(rho_P_back[j]), "S_cn_only": float(S[j]),
                         "pred_kd": "down" if -np.sign(s_ctx[j]) * s_A < 0 else "up", "TRAIN_hit": hit is not None,
                         "observed": None if hit is None else ("down" if hit < 0 else "up"),
                         "concordant": None if hit is None else bool(int(-np.sign(s_ctx[j]) * s_A) == int(hit))})
    ctx_block = {
        "formula": "S_ctx = sign(t_AB) * rho_P * rho_T^2; rho_P = percentile of |S| among all finite trans pairs "
                   "(12,141 x 18,276 matrix, zeros included); rho_T = percentile of B's DepMap HCT116 log2(TPM+1) among 18,276 targets",
        "n_trans_finite_entries": int(N),
        "rho_P_check": [{"target": B, "n_trans_pairs_with_larger_abs_S": int(n_greater[k]),
                         "one_minus_rho_P_exact": float(1 - rho_P_exact[k]),
                         "one_minus_rho_P_from_file": float(1 - rho_P_back[jT[B]])}
                        for k, B in enumerate(q_targets)],
        "n_trans_entries_exactly_zero (G1-failed regulators)": int(n_zero),
        "top10_universe_by_S_ctx": ctx_rows,
    }

    timing["total_s"] = time.time() - T0
    out = {
        "about": "Worked example for docs/explainer/worked_example.md; generated by make_worked_example.py. "
                 "All labels are HCT116 TRAIN Perturb-seq (PIE X-Atlas/Orion, p_adj<0.05).",
        "method": SEL, "params": {k: v for k, v in p.to_dict().items() if k in
                                  ("arm_handling", "cn_floor", "cn_cap", "g1_tau", "g1_w", "a3_rho", "a3_s", "a3_wco", "a3_sigma0")},
        "regulator": {"key": REG, "symbol": sym(REG), "chrom": c, "arm": str(genes.loc[REG, "arm_id"]),
                      "start_bp": int(genes.loc[REG, "start"]), "hct116_split": "TRAIN",
                      "own_class": str(row.own_class.iloc[0]), "cn_tercile": r4(row.cn_tercile.iloc[0]),
                      "expr_in_DepMap_HCT116": float(row.expr_context.iloc[0])},
        "step1_cn_distribution": cn_block,
        "step2_cis_instrument_G1": cis_block,
        "step3_segment": {"rho": p.a3_rho, "segment_index_on_chrom": int(seg), "n_segments_on_chrom": int(labels.max() + 1),
                          "n_cn_genes_on_chrom": int(len(chrom_keys)), "members": seg_rows,
                          "contrast_neighbour": None if k_low is None else members[k_low],
                          "prior_rule": "pi_g proportional to exp(t_cis,g / 2); CN-only genes t_cis = 0",
                          "abf_rule": "log ABF(z) = 0.5*log(1-lam) + 0.5*lam*z^2, lam = 0.05^2/(0.05^2 + 1/df), z = t_gB"},
        "step4_trans_examples": walk,
        "step4_fwl_scatter": fwl_scatter,
        "step5_score": {"rule": "S_AB = t_AB * w_G1(A) * alpha_{A|B}; predicted knockdown sign = -sign(S) * s_A",
                        "w_G1": float(w1), "s_A": s_A,
                        "n_trans_targets": int((trans & np.isfinite(S)).sum()),
                        "abs_t_quantiles_trans": {f"q{int(q*100):02d}": float(np.nanquantile(np.abs(tA[trans]), q)) for q in (0.5, 0.9, 0.99, 1)},
                        "alpha_quantiles_trans": {f"q{int(q*100):02d}": float(np.quantile(aA[trans], q)) for q in (0.01, 0.1, 0.5, 0.9)},
                        "check_vs_committed_scores": check},
        "step6_top10_all_trans_targets": top_all,
        "step6_top10_hct116_universe_targets": top_univ,
        "step6_ranking_and_TRAIN_labels": rank_block,
        "step7_context_layer": ctx_block,
        "timing_s (this script, laptop)": timing,
    }
    tr_pos = {r: k for k, r in enumerate(train_regs.regulator)}
    details = []
    for k in others:
        g = members[k]
        gi = int(np.flatnonzero(regs_all == g)[0]) if g in set(regs_all) else None
        details.append({
            "gene": g, "t_cis": float(cisS[k]), "G1_pass": bool(gates.g1_pass(np.array([cisS[k]]), p.g1_tau)[0]),
            "alpha_median_trans": float(np.median(alphaS[k][trans])),
            "max_abs_t_trans": float(np.nanmax(np.abs(tS[k][trans]))),
            "max_abs_S_in_file": None if gi is None else float(np.nanmax(np.abs(np.asarray(Ssel[gi])))),
            "TRAIN_hits": int((pos.regulator == g).sum()) if g in tr_pos else None,
            "pairs_in_TRAIN_topK": int(ntop[tr_pos[g]]) if g in tr_pos else None,
        })
    out["step3_segment"]["neighbour_details"] = details
    Path(a.out).write_text(json.dumps(r4(out), indent=None, separators=(",", ":")))
    print(f"wrote {a.out} ({Path(a.out).stat().st_size / 1024:.1f} KB) in {timing['total_s']:.1f}s")
    print(json.dumps(r4(check), indent=1))


if __name__ == "__main__":
    main()
