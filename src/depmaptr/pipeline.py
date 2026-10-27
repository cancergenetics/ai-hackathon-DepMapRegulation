"""Variant orchestration: M1 and its pre-registered variants (DESIGN §5.1-§5.6), NC1, B6.

`run(ds, Params)` returns the reduced-form t (regulators x targets), the final score S_AB, cis t, s_A and
diagnostics. All algebra is float64; outputs are cast to float32 by the caller.
"""
from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd

from . import attribution, gates
from .assoc import Residualiser, centred_norm, chunked_pair_partial, diag_partial, lineage_permutation, resid_ss
from .calibrate import lambda_gc
from .covariates import DesignBuilder, arm_cn, build_Z, lineage_halves, lineage_levels, sex_centre, sex_groups
from .io import CHROMS, HOMDEL_MIN_LINES, HOMDEL_RELCN, Dataset, cn_transform, rint

ROOT_SEED = 20261027


@dataclass
class Params:
    arm_handling: str = "loco"          # loco | focal | loao | pcs
    n_pcs: int = 20
    cn_floor: float = 2.0 ** -2
    cn_cap: float = 2.0 ** 3
    rint: bool = False                  # RINT of gene-level x and of y
    dosage: str = "linear"              # linear | hinge (loss slope ranked)
    g1_tau: float | None = None
    g1_w: float = 0.0
    g2_threshold: float | None = None
    g2_w: float = 0.0
    context: str = "HCT116"
    a3_rho: float | None = None
    a3_s: float = 2.0
    a3_wco: float = 0.0
    a3_sigma0: float = 0.05
    target_std: bool = False
    nc1_index: int | None = None        # NC1 residual permutation (SeedSequence(20261027).spawn index)
    sample_set: str = "S2"
    chunk: int = 512
    # --- post-hoc tier D6 (reviews/ideas/M5a_ideas.md §2) ---
    comp_renorm: bool = False           # I5: per group c, y_B -> log2(TPM_B / (1 - share_c) + 1)
    score_base: str = "t"               # t | iv_shrink (I3) | t_rep (I6)
    fuse_gamma: float | None = None     # I1: S = sign(t) * rho_P * rho_T^gamma
    regctx: bool = False                # I2: S <- S * rho_R(A)
    lineage_halves_file: str = "results/split/lineage_halves.tsv"   # I6

    @classmethod
    def from_dict(cls, d: dict) -> "Params":
        known = set(cls.__dataclass_fields__)
        bad = set(d) - known
        if bad:
            raise ValueError(f"unknown params {sorted(bad)}")
        return cls(**d)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class Prepared:
    ds: Dataset
    p: Params
    x_all: np.ndarray            # active CN transform (gene level; RINT'd if p.rint), aligned to cn_keys
    x_lin: np.ndarray            # log2 clip transform (pre-RINT)
    y_all: np.ndarray            # expression (RINT'd if p.rint), aligned to expr_keys
    arms: pd.DataFrame
    Z: np.ndarray
    zlog: dict
    strata: np.ndarray
    cn_idx: dict
    ex_idx: dict
    XT: np.ndarray
    YT: np.ndarray
    DEL: np.ndarray
    y_raw: np.ndarray
    builder: DesignBuilder
    diag: dict = field(default_factory=dict)
    _tpm_tot: np.ndarray | None = None
    _tpm_chr: dict = field(default_factory=dict)

    def group_chrom(self, g: str) -> str:
        """Chromosome excluded by group g (I5). Undefined for the focal variant's single group."""
        if self.p.arm_handling in ("loco", "pcs"):
            return g
        if self.p.arm_handling == "loao":
            return g[:-1]
        raise ValueError("I5 (comp_renorm) is defined per regulator chromosome; not for arm_handling=focal")

    def y_cols(self, g: str, idx) -> np.ndarray:
        """Expression columns `idx` (into expr_keys) as used in group g.

        I5: y -> log2(TPM / (1 - share_c) + 1), TPM = 2^y - 1, share_c = sum_{genes on chr c} TPM / sum_all TPM
        per line, over all (protein-coding) expression columns; RINT (if active) is applied afterwards."""
        idx = np.asarray(idx)
        if not self.p.comp_renorm:
            return self.y_all[:, idx]
        c = self.group_chrom(g)
        if self._tpm_tot is None:
            tpm = np.exp2(np.asarray(self.ds.expr, np.float64)) - 1.0
            self._tpm_tot = tpm.sum(1)
            ch = self.ds.genes.loc[self.ds.expr_keys, "chrom"].astype(str).to_numpy()
            self._tpm_chr = {cc: tpm[:, ch == cc].sum(1) for cc in pd.unique(ch)}
        share = self._tpm_chr.get(c, np.zeros(self.ds.n)) / self._tpm_tot
        tpm = np.exp2(np.asarray(self.ds.expr[:, idx], np.float64)) - 1.0
        y = np.log2(tpm / (1.0 - share)[:, None] + 1.0)
        return rint(y) if self.p.rint else y


def prepare(ds: Dataset, p: Params, Z_df: pd.DataFrame | None = None) -> Prepared:
    """Z_df: optional fixed covariate matrix (rows aligned to ds.panel), used by I6 lineage halves."""
    x_lin = cn_transform(ds.relcn, floor=p.cn_floor, cap=p.cn_cap)
    # Arm CN and burden always use the log2 clip transform (burden thresholds are in log2 units);
    # under RINT only the gene-level regressors and y are rank-normalised (docs/MODEL.md).
    arms = arm_cn(ds, x=x_lin)
    x_all = rint(x_lin) if p.rint else x_lin
    y_all = rint(ds.expr) if p.rint else np.asarray(ds.expr, np.float64)
    if Z_df is None:
        Zdf, zlog = build_Z(ds.panel)
    else:
        Zdf = Z_df
        zlog = {"columns": list(Z_df.columns), "dropped_duplicate": [], "dropped_constant": [],
                "lineage_levels": lineage_levels(ds.panel).value_counts().to_dict()}
    strata = lineage_levels(ds.panel).to_numpy()
    cn_idx = {k: i for i, k in enumerate(ds.cn_keys)}
    ex_idx = {k: i for i, k in enumerate(ds.expr_keys)}
    tc = np.array([cn_idx[k] for k in ds.targets])
    te = np.array([ex_idx[k] for k in ds.targets])
    XT = x_all[:, tc]
    YT = y_all[:, te]
    DEL = (ds.relcn[:, tc] < HOMDEL_RELCN).astype(np.float64)
    ndel = DEL.sum(0)
    DEL[:, ndel < HOMDEL_MIN_LINES] = 0.0
    builder = DesignBuilder(ds, Zdf.to_numpy(), arms, p.arm_handling, x_all=x_all, n_pcs=p.n_pcs)
    prep = Prepared(ds=ds, p=p, x_all=x_all, x_lin=x_lin, y_all=y_all, arms=arms, Z=Zdf.to_numpy(), zlog=zlog,
                    strata=strata, cn_idx=cn_idx, ex_idx=ex_idx, XT=XT, YT=YT, DEL=DEL, y_raw=centred_norm(YT),
                    builder=builder)
    prep.diag.update({"n_targets_with_homdel_indicator": int((ndel >= HOMDEL_MIN_LINES).sum()),
                      "n_arms": int(arms.shape[1]), "arms": list(arms.columns), "Z_columns": zlog["columns"],
                      "Z_dropped_duplicate": zlog["dropped_duplicate"], "Z_dropped_constant": zlog["dropped_constant"],
                      "lineage_levels": {str(k): int(v) for k, v in zlog["lineage_levels"].items()}})
    return prep


class GroupCache:
    """Residualiser + residualised targets/cis basis for one design group (size-1 cache)."""

    def __init__(self, prep: Prepared, builder: DesignBuilder | None = None):
        self.prep, self.builder = prep, builder or prep.builder
        self.key, self.val = None, None
        self.log = {}

    def get(self, g: str, need_targets: bool = True):
        if self.key == (g, need_targets):
            return self.val
        D, names = self.builder.design(g)
        res = Residualiser(D)
        self.log[g] = {"rank": res.rank, "n_cols": res.n_cols, "cond": res.cond}
        val = {"res": res}
        if need_targets:
            pr = self.prep
            if pr.p.comp_renorm:
                te = [pr.ex_idx[k] for k in pr.ds.targets]
                yt = pr.y_cols(g, te)
                val["Y"], val["y_raw"] = res(yt), centred_norm(yt)
            else:
                val["Y"], val["y_raw"] = res(pr.YT), pr.y_raw
            val["W"] = [res(pr.XT), res(pr.DEL)]
        self.key, self.val = (g, need_targets), val
        return val


def nc1_permutation(prep: Prepared, index: int) -> np.ndarray:
    ss = np.random.SeedSequence(ROOT_SEED).spawn(index + 1)[index]
    return lineage_permutation(prep.strata, np.random.default_rng(ss))


def _rows_by_chrom(ds: Dataset, keys) -> dict:
    ch = ds.genes.loc[list(keys), "chrom"].astype(str).to_numpy()
    keys = np.asarray(keys)
    return {c: keys[ch == c] for c in CHROMS if (ch == c).any()}


def cis_t(prep: Prepared, keys, builder: DesignBuilder | None = None) -> np.ndarray:
    """G1 cis t of y_g ~ x_g | D_group(g) (linear x, active transform). 0 for genes without expression."""
    ds = prep.ds
    out = np.zeros(len(keys))
    pos = {k: i for i, k in enumerate(keys)}
    b = builder or prep.builder
    cache = GroupCache(prep, b)
    groups = b.group_of(keys)
    for g in pd.unique(groups):
        ks = np.asarray(keys)[groups == g]
        ks = [k for k in ks if k in prep.ex_idx]
        if not ks:
            continue
        res = cache.get(g, need_targets=False)["res"]
        x = prep.x_all[:, [prep.cn_idx[k] for k in ks]]
        y = prep.y_cols(g, [prep.ex_idx[k] for k in ks])
        _, t, _ = diag_partial(res(x), res(y), df_base=res.df_base, u_raw=centred_norm(x), y_raw=centred_norm(y))
        out[[pos[k] for k in ks]] = np.nan_to_num(t, nan=0.0)
    return out


def group_t(prep: Prepared, gv: dict, keys, perm=None):
    """Reduced-form t, df (len(keys) x T) for regulator rows `keys` of one design group (production path)."""
    p = prep.p
    res, Y, W = gv["res"], gv["Y"], gv["W"]
    x = prep.x_all[:, [prep.cn_idx[k] for k in keys]]
    xg = None
    if p.dosage == "hinge":
        # hinge at copy-neutral: for X regulators that is the sex-group median (male X sits at log2 ~ -1)
        isx = prep.ds.genes.loc[list(keys), "chrom"].astype(str).to_numpy() == "X"
        if isx.any():
            x = x.copy()
            x[:, isx] = sex_centre(x[:, isx], sex_groups(prep.ds.panel))
        xl, xg = np.minimum(x, 0.0), np.maximum(x, 0.0)
        U, V, u_raw = res(xl), [res(xg)], centred_norm(xl)
    else:
        U, V, u_raw = res(x), [], centred_norm(x)
    if perm is not None:
        U = U[perm]
        V = [v[perm] for v in V]
    _, t, df = chunked_pair_partial(U, Y, V, W, df_base=res.df_base, u_raw=u_raw, y_raw=gv["y_raw"], chunk=p.chunk)
    return t, df, (U, V, xg)


def run(ds: Dataset, p: Params, log=print, chroms=None, Z_df: pd.DataFrame | None = None) -> dict:
    """chroms: optional subset of regulator chromosomes to compute (tests); other rows stay NaN."""
    t0 = time.time()
    prep = prepare(ds, p, Z_df=Z_df)
    regs = np.asarray(ds.regulators)
    T = len(ds.targets)
    R = len(regs)
    reg_pos = {k: i for i, k in enumerate(regs)}
    a3 = p.a3_rho is not None
    row_keys = np.asarray(ds.cn_keys) if a3 else regs      # A3 needs t for ALL CN genes in a segment
    by_chrom = _rows_by_chrom(ds, row_keys)
    if chroms is not None:
        by_chrom = {c: v for c, v in by_chrom.items() if c in set(chroms)}
    perm = nc1_permutation(prep, p.nc1_index) if p.nc1_index is not None else None

    t_out = np.full((R, T), np.nan, np.float32)
    alpha_out = np.ones((R, T), np.float32) if a3 else None
    tgain_out = np.full((R, T), np.nan, np.float32) if p.dosage == "hinge" else None
    if p.score_base not in ("t", "iv_shrink", "t_rep"):
        raise ValueError(f"score_base {p.score_base}")
    base_out = np.full((R, T), np.nan, np.float32) if p.score_base == "iv_shrink" else None
    cis_out = np.zeros(R)
    cache = GroupCache(prep)
    loco_builder = DesignBuilder(ds, prep.Z, prep.arms, "loco") if a3 else None
    lam = {}
    ctx = ds.context_expr.get(p.context) if (p.g2_threshold is not None) else None
    if p.g2_threshold is not None and ctx is None:
        raise ValueError(f"context {p.context} has no expression in the dataset")
    n_segments = {}

    for c, keys in by_chrom.items():
        tc0 = time.time()
        m = len(keys)
        t_c = np.full((m, T), np.nan)
        df_c = np.full((m, T), np.nan)
        zco_c = np.zeros((m, T)) if (a3 and p.a3_wco > 0) else None
        dfco_c = np.full((m, T), np.nan) if zco_c is not None else None
        tg_c = np.full((m, T), np.nan) if tgain_out is not None else None
        cis_c = np.zeros(m)
        iv_c = ({"r_cis": np.zeros(m), "sd_yA": np.zeros(m), "sd_B": np.full((m, T), np.nan)}
                if p.score_base == "iv_shrink" else None)
        groups = prep.builder.group_of(keys)
        for g in pd.unique(groups):
            sel = np.flatnonzero(groups == g)
            gv = cache.get(g)
            res, Y, W = gv["res"], gv["Y"], gv["W"]
            t, df, (U, V, xg) = group_t(prep, gv, keys[sel], perm)
            t_c[sel], df_c[sel] = t, df
            if iv_c is not None:   # I3: residual SD of y_B after D_g and W_B (per target)
                iv_c["sd_B"][sel] = (np.sqrt(resid_ss(Y, W)) / np.sqrt(ds.n))[None, :]
            if tg_c is not None:   # gain slope (reported only)
                _, tg, _ = chunked_pair_partial(V[0], Y, [U], W, df_base=res.df_base, u_raw=centred_norm(xg),
                                                y_raw=gv["y_raw"], chunk=p.chunk)
                tg_c[sel] = tg
            # G1 cis t (linear x, unpermuted)
            has_e = np.array([k in prep.ex_idx for k in keys[sel]])
            if has_e.any():
                ks = keys[sel][has_e]
                xe = prep.x_all[:, [prep.cn_idx[k] for k in ks]]
                ye = prep.y_cols(g, [prep.ex_idx[k] for k in ks])
                Xr, Yr = res(xe), res(ye)
                rci, tci, _ = diag_partial(Xr, Yr, df_base=res.df_base, u_raw=centred_norm(xe), y_raw=centred_norm(ye))
                cis_c[sel[has_e]] = np.nan_to_num(tci, nan=0.0)
                if iv_c is not None:
                    iv_c["r_cis"][sel[has_e]] = np.nan_to_num(rci, nan=0.0)
                    iv_c["sd_yA"][sel[has_e]] = np.linalg.norm(Yr, axis=0) / np.sqrt(ds.n)
                if zco_c is not None:
                    rco, _, dfco = chunked_pair_partial(Yr, Y, [Xr], W, df_base=res.df_base,
                                                        u_raw=centred_norm(ye), y_raw=gv["y_raw"], chunk=p.chunk)
                    with np.errstate(invalid="ignore", divide="ignore"):
                        z = np.arctanh(np.clip(rco, -1 + 1e-15, 1 - 1e-15)) * np.sqrt(dfco - 1)
                    zco_c[sel[has_e]] = np.nan_to_num(z, nan=0.0)
                    dfco_c[sel[has_e]] = dfco
        if dfco_c is not None:
            # CN-only genes (no expression column) have no co-expression evidence: z_co = 0, with the df the
            # co-expression model would have (one more regressor, x_g). Never NaN (PI decision D5).
            dfco_c = np.where(np.isfinite(dfco_c), dfco_c, df_c - 1)
        if a3:
            Dl, _ = loco_builder.design(c)
            xres = Residualiser(Dl)(prep.x_all[:, [prep.cn_idx[k] for k in keys]])
            labels = attribution.greedy_segments(xres, p.a3_rho)
            n_segments[c] = int(labels.max() + 1)
            wg2 = gates.g2_weight(keys, ctx, p.g2_threshold, p.g2_w) if ctx is not None else None
            alpha = attribution.segment_alpha(labels, t_c, df_c, cis_c, p.a3_s, w_g2=wg2, zco=zco_c,
                                              w_co=p.a3_wco, df_co=dfco_c if dfco_c is not None else None,
                                              sigma0=p.a3_sigma0)
        rr = np.array([reg_pos.get(k, -1) for k in keys])
        keep = rr >= 0
        t_out[rr[keep]] = t_c[keep]
        if iv_c is not None:
            base_out[rr[keep]] = iv_shrink(t_c[keep], df_c[keep], iv_c["sd_B"][keep], iv_c["r_cis"][keep],
                                           iv_c["sd_yA"][keep])
        cis_out[rr[keep]] = cis_c[keep]
        if a3:
            alpha_out[rr[keep]] = alpha[keep]
        if tg_c is not None:
            tgain_out[rr[keep]] = tg_c[keep]
        # lambda_GC on trans pairs of this regulator chromosome (subsample)
        tchrom = ds.genes.loc[ds.targets, "chrom"].astype(str).to_numpy()
        tt = t_c[keep][:, tchrom != c]
        dd = df_c[keep][:, tchrom != c]
        if tt.size:
            rng = np.random.default_rng(0)
            idx = rng.choice(tt.size, size=min(tt.size, 200_000), replace=False)
            lam[c] = lambda_gc(tt.ravel()[idx], dd.ravel()[idx])
        log(f"  chr{c}: {keep.sum()} regulators ({m} rows), {time.time() - tc0:.1f}s, lambda_GC={lam.get(c, np.nan):.3f}")

    # A = B entries are NaN in every run (they are finite under NC1 otherwise)
    tpos = {k: j for j, k in enumerate(ds.targets)}
    same = np.array([(i, tpos[k]) for i, k in enumerate(regs) if k in tpos], dtype=np.int64).reshape(-1, 2)
    if len(same):
        t_out[same[:, 0], same[:, 1]] = np.nan
    extras_pre = {}
    if p.score_base == "t_rep":
        if p.nc1_index is not None:
            raise ValueError("t_rep (I6) is not defined together with NC1")
        base_out, extras_pre = lineage_replicated_t(ds, p, Z_df, log=log, chroms=chroms)
        base_out[np.isnan(t_out)] = np.nan
    if base_out is not None and len(same):
        base_out[same[:, 0], same[:, 1]] = np.nan
    score, extras = compose_score(ds, p, t_out, cis_out, alpha_out, base=base_out)
    arr_extras = {k[len("_arr_"):]: extras.pop(k) for k in [k for k in extras if k.startswith("_arr_")]}
    extras.update(extras_pre)
    new_nan = np.isnan(score) & ~np.isnan(t_out)
    assert not new_nan.any(), f"score has {int(new_nan.sum())} NaNs where t is finite"
    diag = dict(prep.diag)
    diag.update({"groups": cache.log, "lambda_gc_by_chrom": lam, "n_segments_by_chrom": n_segments,
                 "wall_engine_s": time.time() - t0, "nc1_index": p.nc1_index, **extras})
    out = {"t": t_out, "score": score, "cis_t": cis_out.astype(np.float32), "s_cis": gates.s_cis(cis_out),
           "diagnostics": diag}
    if alpha_out is not None:
        out["alpha"] = alpha_out
    if tgain_out is not None:
        out["t_gain"] = tgain_out
    out.update(arr_extras)
    return out


def target_standardise(ds: Dataset, t: np.ndarray) -> np.ndarray:
    """z-score each target's t across regulators, excluding same-chromosome (and NaN) entries."""
    rch = ds.genes.loc[ds.regulators, "chrom"].astype(str).to_numpy()
    tch = ds.genes.loc[ds.targets, "chrom"].astype(str).to_numpy()
    out = np.empty_like(t)
    for c in pd.unique(tch):
        cols = np.flatnonzero(tch == c)
        blk = t[:, cols].astype(np.float64)
        trans = blk[rch != c]
        mu = np.nanmean(trans, axis=0)
        sd = np.nanstd(trans, axis=0, ddof=1)
        out[:, cols] = ((blk - mu) / sd).astype(t.dtype)
    return out


def compose_score(ds: Dataset, p: Params, t: np.ndarray, cis: np.ndarray, alpha: np.ndarray | None,
                  base: np.ndarray | None = None):
    """S_inc = b_AB * w_G1(A) * w_G2(A) * alpha_{A|B} (DESIGN §5.6 as amended by PI decision D5: w_G2 is always
    applied as a multiplier, including under A3 where it also enters the prior pi_g).

    b = t (default), the shrunk IV effect (I3) or the lineage-replicated t (I6); target-standardised if active.
    Post-hoc D6 layers, in this order: I1 fusion S = sign(t) * rho_P(|S_inc|) * rho_T^gamma; I2 S <- S * rho_R(A).
    """
    extras = {}
    base = (t if base is None else base)
    base = target_standardise(ds, base) if p.target_std else base.copy()
    w = np.ones(len(ds.regulators))
    if p.g1_tau is not None:
        w1 = gates.g1_weight(cis, p.g1_tau, p.g1_w)
        extras["n_g1_pass"] = int((w1 == 1).sum())
        w = w * w1
    if p.g2_threshold is not None:
        ctx = ds.context_expr[p.context]
        expressed = gates.context_expressed(ds.regulators, ctx, p.g2_threshold)
        extras["n_g2_expressed"] = int(expressed.sum())
        # PI decision D5: G2 is ALWAYS a multiplier (also under A3, where it additionally enters pi)
        w = w * np.where(expressed, 1.0, p.g2_w)
    score = base * w[:, None].astype(np.float32)
    if alpha is not None:
        score = score * alpha
    score = score.astype(np.float32)
    if p.fuse_gamma is not None:
        score, ex = context_fusion(ds, p, t, score)
        extras.update(ex)
    if p.regctx:
        rho_r = regulator_context_percentile(ds, p.context)
        extras["regctx_n_zero"] = int((rho_r == 0).sum())
        score = (score * rho_r[:, None].astype(np.float32)).astype(np.float32)
        if "_arr_pair_score" in extras:
            extras["_arr_pair_score"] = (extras["_arr_pair_score"] * rho_r[:, None].astype(np.float32)).astype(np.float32)
    return score, extras


# ----------------------------------------------------------------------------- D6 post-hoc layers
def iv_shrink(t, df, sd_B, r_cis, sd_yA):
    """I3: beta_IV = r_AB * sd_B / (r_cis * sd_yA), shrunk by max(0, 1 - 4/t^2) (positive-part James-Stein).

    r_AB = t / sqrt(df + t^2) (the engine's partial correlation). Regulators without a positive cis correlation
    (incl. CN-only, r_cis = 0) have no instrument: beta = 0 (docs/MODEL.md §10)."""
    t = np.asarray(t, np.float64)
    r = t / np.sqrt(df + t * t)
    ok = (np.asarray(r_cis) > 0) & (np.asarray(sd_yA) > 0)
    denom = np.where(ok, r_cis * sd_yA, 1.0)[:, None]
    beta = np.where(ok[:, None], r * sd_B / denom, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        shrink = np.maximum(0.0, 1.0 - 4.0 / (t * t))
    return (beta * shrink).astype(np.float32)


def trans_percentile(ds: Dataset, v: np.ndarray) -> np.ndarray:
    """Average-rank percentile (rank / N, ties averaged) of v among all finite TRANS entries of the matrix.
    Non-trans entries get 0; NaN stays NaN. One argsort over the trans entries (O(N log N))."""
    rch = ds.genes.loc[ds.regulators, "chrom"].astype(str).to_numpy()
    tch = ds.genes.loc[ds.targets, "chrom"].astype(str).to_numpy()
    mask = (rch[:, None] != tch[None, :]) & np.isfinite(v)
    vals = v[mask]
    N = len(vals)
    order = np.argsort(vals, kind="stable")
    sv = vals[order]
    del vals
    starts = np.flatnonzero(np.r_[True, sv[1:] != sv[:-1]])
    ends = np.r_[starts[1:], N]
    del sv
    pct_sorted = np.repeat(((starts + ends + 1) / 2.0) / N, ends - starts)     # average 1-based rank / N
    pct = np.empty(N, np.float32)
    pct[order] = pct_sorted
    del order, pct_sorted
    out = np.zeros(v.shape, np.float32)
    out[mask] = pct
    out[~np.isfinite(v)] = np.nan
    return out


def target_context_percentile(ds: Dataset, context: str) -> np.ndarray:
    """rho_T(B): average-rank percentile of context log2(TPM+1) among ALL measured targets (never the PIE universe)."""
    from scipy.stats import rankdata
    v = ds.context_expr[context].reindex(list(ds.targets)).to_numpy(np.float64)
    assert np.isfinite(v).all(), "context expression missing for some targets"
    return rankdata(v) / len(v)


def regulator_context_percentile(ds: Dataset, context: str) -> np.ndarray:
    """rho_R(A): average-rank percentile of context log2(TPM+1) among the tool regulators that have an
    expression column; CN-only regulators get 0."""
    from scipy.stats import rankdata
    v = ds.context_expr[context].reindex(list(ds.regulators)).to_numpy(np.float64)
    ok = np.isfinite(v)
    out = np.zeros(len(v))
    out[ok] = rankdata(v[ok]) / ok.sum()
    return out


def context_fusion(ds: Dataset, p: Params, t: np.ndarray, s_inc: np.ndarray):
    """I1: S = sign(t_AB) * rho_P(A,B) * rho_T(B)^gamma; rho_P = trans percentile of |S_inc| (0 where S_inc = 0)."""
    a = np.abs(s_inc)
    rho_p = trans_percentile(ds, a)
    rho_p[a == 0] = 0.0
    rho_t = target_context_percentile(ds, p.context) ** p.fuse_gamma
    with np.errstate(invalid="ignore"):
        pair = (np.sign(t) * rho_p).astype(np.float32)
    pair[np.isnan(s_inc)] = np.nan
    S = pair * rho_t[None, :].astype(np.float32)
    # optional INTERFACES fields for the exact C1 control: score == pair_score * target_weight
    return S.astype(np.float32), {"fusion_gamma": p.fuse_gamma, "fusion_n_zero_rho_p": int((a == 0).sum()),
                                  "_arr_pair_score": pair, "_arr_target_weight": rho_t.astype(np.float32)}


def load_lineage_halves(ds: Dataset, p: Params) -> pd.Series:
    """I6 split: committed file for real data (asserted to cover exactly the panel's lineage levels);
    recomputed with the same seeded rule for SIM."""
    from pathlib import Path
    from .io import ROOT
    levels = lineage_levels(ds.panel)
    f = Path(p.lineage_halves_file)
    f = f if f.is_absolute() else ROOT / f
    if ds.sample_set == "SIM":
        split = lineage_halves(ds.panel)
    else:
        tab = pd.read_csv(f, sep="\t")
        split = tab.set_index("lineage_level")["half"]
        if set(split.index) != set(levels.unique()):
            raise AssertionError(f"{f} does not match the panel's lineage levels")
    return levels.map(split)


def lineage_replicated_t(ds: Dataset, p: Params, Z_df=None, log=print, chroms=None):
    """I6: t_rep = sign(t1) min(|t1|, |t2|) 1[sign t1 = sign t2], M1 fitted separately in each lineage half
    (full-sample Z restricted to the half, constant/duplicated columns dropped). A half-NaN gives t_rep = 0."""
    half = load_lineage_halves(ds, p).to_numpy()
    Zfull = build_Z(ds.panel)[0] if Z_df is None else Z_df
    ph = Params.from_dict({**p.to_dict(), "score_base": "t", "fuse_gamma": None, "regctx": False, "a3_rho": None,
                           "g1_tau": None, "g2_threshold": None, "target_std": False})
    ts, info = [], {}
    for h in (1, 2):
        m = half == h
        dsh = ds.subset_rows(m, f"{ds.sample_set}_H{h}")
        Zh = Zfull.loc[m].reset_index(drop=True)
        Zh = Zh.loc[:, Zh.nunique() > 1]
        Zh = Zh.loc[:, ~Zh.T.duplicated().to_numpy()]
        log(f"  I6 lineage half {h}: n={int(m.sum())}, Z columns {Zh.shape[1]}")
        ts.append(run(dsh, ph, log=lambda *a, **k: None, chroms=chroms, Z_df=Zh)["t"])
        info[f"half{h}_n"] = int(m.sum())
        info[f"half{h}_Z_columns"] = int(Zh.shape[1])
    t1, t2 = ts
    with np.errstate(invalid="ignore"):
        rep = np.sign(t1) * np.minimum(np.abs(t1), np.abs(t2)) * (np.sign(t1) == np.sign(t2))
    rep = np.where(np.isnan(rep), 0.0, rep).astype(np.float32)
    info["t_rep_frac_sign_agree"] = float(np.nanmean(np.sign(t1) == np.sign(t2)))
    return rep, info


# ----------------------------------------------------------------------------- B6 arm-only
def run_arm_only(ds: Dataset, p: Params, log=print) -> dict:
    """B6: t of (arm CN of A's arm) on y_B given D_c (LOCO) and W_B; broadcast to regulators on that arm."""
    prep = prepare(ds, p)
    regs = np.asarray(ds.regulators)
    T = len(ds.targets)
    t_out = np.full((len(regs), T), np.nan, np.float32)
    reg_arm = ds.genes.loc[regs, "arm_id"].astype(str).to_numpy()
    reg_chr = ds.genes.loc[regs, "chrom"].astype(str).to_numpy()
    cache = GroupCache(prep)
    for c in pd.unique(reg_chr):
        arms_c = [a for a in prep.arms.columns if a[:-1] == c]
        if not arms_c:
            continue
        gv = cache.get(c)
        res = gv["res"]
        A = prep.arms[arms_c].to_numpy()
        _, t, _ = chunked_pair_partial(res(A), gv["Y"], [], gv["W"], df_base=res.df_base, u_raw=centred_norm(A),
                                       y_raw=gv["y_raw"])
        for j, a in enumerate(arms_c):
            t_out[reg_arm == a] = t[j]
        log(f"  chr{c}: arms {arms_c}")
    cis = cis_t(prep, regs)
    return {"t": t_out, "score": t_out.copy(), "cis_t": cis.astype(np.float32), "s_cis": gates.s_cis(cis),
            "diagnostics": {**prep.diag, "groups": cache.log}}
