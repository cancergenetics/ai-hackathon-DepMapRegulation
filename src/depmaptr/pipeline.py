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
from .assoc import Residualiser, centred_norm, chunked_pair_partial, diag_partial, lineage_permutation
from .calibrate import lambda_gc
from .covariates import DesignBuilder, arm_cn, build_Z, lineage_levels, sex_centre, sex_groups
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


def prepare(ds: Dataset, p: Params) -> Prepared:
    x_lin = cn_transform(ds.relcn, floor=p.cn_floor, cap=p.cn_cap)
    # Arm CN and burden always use the log2 clip transform (burden thresholds are in log2 units);
    # under RINT only the gene-level regressors and y are rank-normalised (docs/MODEL.md).
    arms = arm_cn(ds, x=x_lin)
    x_all = rint(x_lin) if p.rint else x_lin
    y_all = rint(ds.expr) if p.rint else np.asarray(ds.expr, np.float64)
    Zdf, zlog = build_Z(ds.panel)
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
            val["Y"] = res(pr.YT)
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
        y = prep.y_all[:, [prep.ex_idx[k] for k in ks]]
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
    _, t, df = chunked_pair_partial(U, Y, V, W, df_base=res.df_base, u_raw=u_raw, y_raw=prep.y_raw, chunk=p.chunk)
    return t, df, (U, V, xg)


def run(ds: Dataset, p: Params, log=print, chroms=None) -> dict:
    """chroms: optional subset of regulator chromosomes to compute (tests); other rows stay NaN."""
    t0 = time.time()
    prep = prepare(ds, p)
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
        groups = prep.builder.group_of(keys)
        for g in pd.unique(groups):
            sel = np.flatnonzero(groups == g)
            gv = cache.get(g)
            res, Y, W = gv["res"], gv["Y"], gv["W"]
            t, df, (U, V, xg) = group_t(prep, gv, keys[sel], perm)
            t_c[sel], df_c[sel] = t, df
            if tg_c is not None:   # gain slope (reported only)
                _, tg, _ = chunked_pair_partial(V[0], Y, [U], W, df_base=res.df_base, u_raw=centred_norm(xg),
                                                y_raw=prep.y_raw, chunk=p.chunk)
                tg_c[sel] = tg
            # G1 cis t (linear x, unpermuted)
            has_e = np.array([k in prep.ex_idx for k in keys[sel]])
            if has_e.any():
                ks = keys[sel][has_e]
                xe = prep.x_all[:, [prep.cn_idx[k] for k in ks]]
                ye = prep.y_all[:, [prep.ex_idx[k] for k in ks]]
                Xr, Yr = res(xe), res(ye)
                _, tci, _ = diag_partial(Xr, Yr, df_base=res.df_base, u_raw=centred_norm(xe), y_raw=centred_norm(ye))
                cis_c[sel[has_e]] = np.nan_to_num(tci, nan=0.0)
                if zco_c is not None:
                    rco, _, dfco = chunked_pair_partial(Yr, Y, [Xr], W, df_base=res.df_base,
                                                        u_raw=centred_norm(ye), y_raw=prep.y_raw, chunk=p.chunk)
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
    score, extras = compose_score(ds, p, t_out, cis_out, alpha_out)
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


def compose_score(ds: Dataset, p: Params, t: np.ndarray, cis: np.ndarray, alpha: np.ndarray | None):
    """S_AB = t_AB * w_G1(A) * w_G2(A) * alpha_{A|B} (DESIGN §5.6 as amended by PI decision D5: w_G2 is always
    applied as a multiplier, including under A3 where it also enters the prior pi_g)."""
    extras = {}
    base = target_standardise(ds, t) if p.target_std else t.copy()
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
    return score.astype(np.float32), extras


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
                                       y_raw=prep.y_raw)
        for j, a in enumerate(arms_c):
            t_out[reg_arm == a] = t[j]
        log(f"  chr{c}: arms {arms_c}")
    cis = cis_t(prep, regs)
    return {"t": t_out, "score": t_out.copy(), "cis_t": cis.astype(np.float32), "s_cis": gates.s_cis(cis),
            "diagnostics": {**prep.diag, "groups": cache.log}}
