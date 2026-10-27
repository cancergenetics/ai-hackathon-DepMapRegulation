"""Shared covariates Z, arm CN, LOCO burden and per-group designs D (DESIGN §4, §5.1, §5.5 tier 1).

Z (no intercept; the intercept is added in `design`):
- lineage: OncotreeLineage; the heterogeneous lineages are split by primary disease where that level has
  >= 15 lines (remainder stays in the parent lineage); lineages with < 15 lines (and missing) -> "other";
- ploidy, WGD: lineage-median imputation + missingness indicator;
- MSI: MSIScore > 20 flag + missingness indicator; TP53 altered + missingness indicator;
- sex: Female / Unknown-or-missing indicators (reference Male);
- media base: levels >= 15 one-hot, rest pooled into "other"; serum-free flag + missingness indicator;
- growth pattern: Adherent reference, other levels (incl. "Unknown" as a missingness level) one-hot;
- RNA technical: Stranded. Platform is dropped (PI decision D3: 93% missing).
Exactly duplicated and constant columns are removed (e.g. wgd_missing == ploidy_missing).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .io import CHROMS, ARM_MIN_GENES, Dataset, cn_transform

HETEROGENEOUS = ("Lymphoid", "Myeloid", "CNS/Brain", "Esophagus/Stomach", "Soft Tissue", "Bone")
MIN_LEVEL = 15
MSI_THRESHOLD = 20.0
BURDEN_THRESHOLD = 0.3
AUTOSOMES = [str(i) for i in range(1, 23)]


def _bool(s: pd.Series) -> pd.Series:
    """Parse True/False/NaN columns that may arrive as object strings."""
    m = {True: 1.0, False: 0.0, "True": 1.0, "False": 0.0, "true": 1.0, "false": 0.0, 1: 1.0, 0: 0.0}
    return s.map(lambda v: m.get(v, np.nan) if not (isinstance(v, float) and np.isnan(v)) else np.nan).astype(float)


def lineage_levels(panel: pd.DataFrame) -> pd.Series:
    lin = panel["lineage"].astype(object).where(panel["lineage"].notna(), None)
    dis = panel["primary_disease"].astype(object)
    out = lin.copy()
    for L in HETEROGENEOUS:
        m = lin == L
        vc = dis[m].value_counts()
        for d, k in vc.items():
            if k >= MIN_LEVEL:
                out[m & (dis == d)] = f"{L}::{d}"
    out = out.where(out.notna(), "other")
    vc = out.value_counts()
    small = set(vc.index[vc < MIN_LEVEL])
    # Remainders of split lineages stay as the parent lineage even if < 15 (they ARE the parent level);
    # whole lineages with < 15 lines go to "other".
    is_parent_remainder = out.isin(HETEROGENEOUS) & lin.isin(HETEROGENEOUS)
    out = out.where(~(out.isin(small) & ~is_parent_remainder), "other")
    return out.rename("lineage_level")


def _onehot(s: pd.Series, prefix: str, ref: str | None = None) -> pd.DataFrame:
    s = s.astype(str)
    levels = s.value_counts()
    ref = ref if ref is not None else levels.index[0]
    cols = {f"{prefix}[{lv}]": (s == lv).astype(float).to_numpy() for lv in sorted(levels.index) if lv != ref}
    return pd.DataFrame(cols, index=s.index)


def build_Z(panel: pd.DataFrame, drop_platform: bool = True) -> tuple[pd.DataFrame, dict]:
    """Return Z (n x p, float64, no intercept) and a log dict."""
    p = panel.reset_index(drop=True)
    log = {}
    lin = lineage_levels(p)
    log["lineage_levels"] = lin.value_counts().to_dict()
    parts = [_onehot(lin, "lineage")]

    def impute(col, miss_col):
        v = p[col].astype(float).copy()
        miss = v.isna() if miss_col is None else p[miss_col].astype(bool)
        med_by = v.groupby(lin).median()
        fill = lin.map(med_by).astype(float).fillna(v.median())
        v = v.where(~miss, fill)
        return v, miss.astype(float)

    ploidy, ploidy_miss = impute("ploidy", "ploidy_missing")
    wgd, wgd_miss = impute("wgd", "wgd_missing")
    msi_miss = p["msi_score"].isna()
    msi = (p["msi_score"].astype(float) > MSI_THRESHOLD).astype(float).where(~msi_miss, 0.0)
    if "msi_high" in p:
        mh = _bool(p["msi_high"])
        ok = mh.notna()
        assert (mh[ok] == msi[ok]).all(), "msi_high disagrees with MSIScore > 20"
    tp = _bool(p["tp53_altered"])
    tp_miss = tp.isna()
    sex = p["sex"].astype(object)
    sf = _bool(p["serum_free"])
    num = pd.DataFrame({
        "ploidy": ploidy, "ploidy_missing": ploidy_miss, "wgd": wgd, "wgd_missing": wgd_miss,
        "msi_high": msi, "msi_missing": msi_miss.astype(float),
        "tp53_altered": tp.fillna(0.0), "tp53_missing": tp_miss.astype(float),
        "sex[Female]": (sex == "Female").astype(float),
        "sex[Unknown]": (~sex.isin(["Female", "Male"])).astype(float),
        "serum_free": sf.fillna(0.0), "serum_free_missing": sf.isna().astype(float),
        "stranded": _bool(p["stranded"]).fillna(0.0),
    })
    parts.append(num)
    media = p["media_base"].astype(object).where(p["media_base"].notna(), "other")
    vc = media.value_counts()
    media = media.where(media.isin(vc.index[vc >= MIN_LEVEL]) & (media != "other"), "other")
    log["media_levels"] = media.value_counts().to_dict()
    parts.append(_onehot(media, "media"))
    gp = p["growth_pattern"].astype(object).where(p["growth_pattern"].notna(), "Unknown")
    parts.append(_onehot(gp, "growth", ref="Adherent" if (gp == "Adherent").any() else None))
    if not drop_platform and "platform" in p:
        plat = p["platform"].astype(object).where(p["platform"].notna(), "missing")
        vc = plat.value_counts()
        plat = plat.where(plat.isin(vc.index[vc >= MIN_LEVEL]), "missing")
        parts.append(_onehot(plat, "platform", ref="missing"))
    Z = pd.concat(parts, axis=1)
    # drop constant and exactly duplicated columns (logged)
    const = [c for c in Z.columns if Z[c].nunique() <= 1]
    Z = Z.drop(columns=const)
    dup = Z.T.duplicated()
    log["dropped_constant"] = const
    log["dropped_duplicate"] = list(Z.columns[dup.to_numpy()])
    Z = Z.loc[:, ~dup.to_numpy()]
    log["columns"] = list(Z.columns)
    return Z.astype(np.float64), log


# ----------------------------------------------------------------------------- arms
def sex_groups(panel: pd.DataFrame) -> np.ndarray:
    s = panel["sex"].astype(object)
    return np.where(s.isin(["Male", "Female"]), s, "Unknown").astype(str)


def sex_centre(M: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Subtract each column's median within sex group (DESIGN §3, §4)."""
    M = np.array(M, np.float64, copy=True)
    for g in np.unique(groups):
        m = groups == g
        M[m] -= np.median(M[m], axis=0)
    return M


def arm_cn(ds: Dataset, x: np.ndarray | None = None, floor: float | None = None) -> pd.DataFrame:
    """arm_k = median of x over CN genes on arm k (arms with >= 20 genes); Xp/Xq sex-centred.

    `x` is the active-transform CN matrix aligned to ds.cn_keys (default log2 clip(2^-2, 2^3)).
    """
    if x is None:
        x = cn_transform(ds.relcn) if floor is None else cn_transform(ds.relcn, floor=floor)
    arm_id = ds.genes.loc[ds.cn_keys, "arm_id"].to_numpy()
    chrom = ds.genes.loc[ds.cn_keys, "chrom"].astype(str).to_numpy()
    cols = {}
    for c in CHROMS:
        for a in ("p", "q"):
            aid = f"{c}{a}"
            m = arm_id == aid
            if m.sum() >= ARM_MIN_GENES:
                cols[aid] = np.median(x[:, m], axis=1)
    A = pd.DataFrame(cols)
    xs = [c for c in A.columns if c.startswith("X")]
    if xs:
        A[xs] = sex_centre(A[xs].to_numpy(), sex_groups(ds.panel))
    return A


def arm_chrom(arm_id: str) -> str:
    return arm_id[:-1]


def burden(arms: pd.DataFrame, exclude: set[str]) -> np.ndarray:
    """LOCO genome burden over AUTOSOMAL arms not in `exclude`: [fraction |arm| > 0.3, mean |arm|]."""
    use = [a for a in arms.columns if arm_chrom(a) in AUTOSOMES and a not in exclude]
    V = np.abs(arms[use].to_numpy())
    return np.column_stack([(V > BURDEN_THRESHOLD).mean(1), V.mean(1)])


# ----------------------------------------------------------------------------- groups and designs
class DesignBuilder:
    """Builds D_g = [1, Z, burden, arm block] for each regulator group of an arm-handling mode.

    Modes (DESIGN §5.5 tier 1):
      loco  : group = chromosome c; arms not on c; burden over autosomal arms not on c.
      focal : one group; all arms (incl. A's own); burden over all autosomal arms.
      loao  : group = arm a; all arms except a; burden over autosomal arms except a.
      pcs   : group = chromosome c; burden LOCO(c); top-k PCs of gene CN on chromosomes != c
              (X genes sex-centred) instead of arms.
    """

    def __init__(self, ds: Dataset, Z: np.ndarray, arms: pd.DataFrame, mode: str = "loco",
                 x_all: np.ndarray | None = None, n_pcs: int = 20, z_names: list[str] | None = None):
        self.ds, self.Z, self.arms, self.mode, self.n_pcs = ds, Z, arms, mode, n_pcs
        self.z_names = list(z_names) if z_names is not None else [f"Z{i}" for i in range(Z.shape[1])]
        assert len(self.z_names) == Z.shape[1]
        if mode == "pcs":
            assert x_all is not None
            chrom = ds.genes.loc[ds.cn_keys, "chrom"].astype(str).to_numpy()
            xc = np.array(x_all, np.float64, copy=True)
            isx = chrom == "X"
            if isx.any():
                xc[:, isx] = sex_centre(xc[:, isx], sex_groups(ds.panel))
            xc -= xc.mean(0)
            self._xc, self._chrom = xc, chrom
            self._gram = xc @ xc.T

    def group_of(self, gene_keys) -> np.ndarray:
        g = self.ds.genes.loc[list(gene_keys)]
        if self.mode in ("loco", "pcs"):
            return g["chrom"].astype(str).to_numpy()
        if self.mode == "focal":
            return np.full(len(g), "all", dtype=object)
        if self.mode == "loao":
            return g["arm_id"].astype(str).to_numpy()
        raise ValueError(self.mode)

    def design(self, group: str) -> tuple[np.ndarray, list[str]]:
        n = self.ds.n
        arms = self.arms
        if self.mode in ("loco", "pcs"):
            excl = {a for a in arms.columns if arm_chrom(a) == group}
        elif self.mode == "focal":
            excl = set()
        else:
            excl = {group}
        B = burden(arms, excl)
        if self.mode == "pcs":
            m = self._chrom == group
            G = self._gram - self._xc[:, m] @ self._xc[:, m].T
            w, U = np.linalg.eigh(G)
            idx = np.argsort(w)[::-1][: self.n_pcs]
            blk = U[:, idx] * np.sqrt(np.maximum(w[idx], 0))
            names = [f"cnpc{i+1}" for i in range(blk.shape[1])]
        else:
            use = [a for a in arms.columns if a not in excl]
            blk = arms[use].to_numpy()
            names = [f"arm[{a}]" for a in use]
        D = np.column_stack([np.ones(n), self.Z, B, blk])
        return D, ["intercept"] + [f"Z:{z}" for z in self.z_names] + ["burden_frac", "burden_mean"] + names


def lineage_halves(panel: pd.DataFrame, seed: int = 20261027) -> pd.Series:
    """I6 (D6): split the Z lineage levels into halves 1/2 by seeded greedy balance on n.

    Levels are taken in decreasing n (ties in a seeded random order) and each goes to the half with the
    smaller current total (exact ties: seeded coin). Returns a Series lineage_level -> half."""
    vc = lineage_levels(panel).value_counts()
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    key = rng.random(len(vc))
    order = sorted(range(len(vc)), key=lambda i: (-int(vc.iloc[i]), key[i]))
    tot = {1: 0, 2: 0}
    out = {}
    for i in order:
        if tot[1] == tot[2]:
            h = int(rng.integers(1, 3))
        else:
            h = 1 if tot[1] < tot[2] else 2
        out[vc.index[i]] = h
        tot[h] += int(vc.iloc[i])
    return pd.Series(out, name="half")
