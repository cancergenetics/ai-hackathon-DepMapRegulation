"""Benchmark harness: universe, labels, held-out split, sparse metrics, resampling, split guard.

Spec: DESIGN.md v3.1 §6 (benchmark definitions and metrics), §7 (baselines use this harness),
§8 (rules: tuning log, single TEST run, invariants); contracts in docs/INTERFACES.md; implementation
choices in docs/BENCHMARK.md.

Benchmark labels are used ONLY for evaluation. Nothing here feeds back into a score.

Layout
------
A. constants and small utilities (hashing, ranks, odds ratios, CMH)
B. building the benchmark tables (pure functions on DataFrames; IO lives in
   scripts/benchmark/build_benchmark.py)
C. held-out split (arm triplets, seed 20261027) and TRAIN folds
D. loading the benchmark and score files, with the TEST guard
E. the evaluation engine (`Evaluation`): top-K selection, pooled OR, permutation p, jackknife,
   two-way bootstrap, CMH, AUROC, sign metrics, call-set OR, hub exclusion
F. `evaluate()` / `compare()` with the automatic tuning log
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

# =====================================================================================
# A. constants and utilities
# =====================================================================================

ROOT = Path(os.environ.get("DEPMAPTR_ROOT", Path(__file__).resolve().parents[2]))
SEED = 20261027                       # DESIGN §8.9; the split uses this seed only
LINES = ("hct116", "hek293t")
PRIMARY_LINE = "hct116"
PRIMARY_K = 1e-3                      # top 0.1 % of universe pairs
K_FRACS = (1e-3, 1e-4, 1e-2)          # 0.1 % (primary), 0.01 %, 1 %
PADJ_POS = 0.05
PADJ_STRICT = 0.01
UEFFECT_MIN = 0.2                     # |Ueffect| >= 0.2 (PIE Ueffect is rank-biserial in [-1,1]; D4)
CONSISTENT_REF_MEAN = 0.05            # PIE control mean >= 0.05 -> tested for every perturbation
BIN_BP = 5_000_000
CONTEXT_MODEL = {"hct116": "ACH-000971", "hek293t": "ACH-001085"}   # DepMap HCT 116 / HEK 293
CHROMS = [str(i) for i in range(1, 23)] + ["X"]
POS_FLAGS = ("pos_primary", "pos_ueffect", "pos_p01")
N_PERM = 1000
N_BOOT = 1000
# numpy.random.SeedSequence(SEED).spawn(N_SPAWN)[index] -> each stochastic component (DESIGN §8.9)
N_SPAWN = 8
SPAWN = {"perm_regulator": 0, "perm_degree": 1, "perm_sign": 2, "bootstrap": 3}
SPLITS = ("TRAIN", "TEST")


def spawn_rng(index: int, seed: int = SEED) -> np.random.Generator:
    """Generator for spawn index `index` of SeedSequence(seed) (DESIGN §8.9)."""
    return np.random.default_rng(np.random.SeedSequence(seed).spawn(N_SPAWN)[index])


def chrom_order(c) -> int:
    c = str(c).removeprefix("chr")
    if c == "X":
        return 23
    if c == "Y":
        return 24
    try:
        return int(c)
    except ValueError:
        return 99


def arm_sort_key(arm_id: str) -> tuple:
    """Genomic order of an arm id like '1p', '17q', 'Xp'."""
    return (chrom_order(arm_id[:-1]), 0 if arm_id[-1] == "p" else 1)


def key_hash(keys: Iterable[str], seed: int = SEED) -> np.ndarray:
    """Stable seeded 64-bit hash of gene keys (blake2b)."""
    out = [int.from_bytes(hashlib.blake2b(f"{seed}|{k}".encode(), digest_size=8).digest(), "little")
           for k in keys]
    return np.array(out, dtype=np.uint64)


def _mix64(x: np.ndarray) -> np.ndarray:
    """splitmix64 finaliser on uint64 arrays (wrapping arithmetic)."""
    x = x.astype(np.uint64, copy=True)
    with np.errstate(over="ignore"):
        x ^= x >> np.uint64(30)
        x *= np.uint64(0xBF58476D1CE4E5B9)
        x ^= x >> np.uint64(27)
        x *= np.uint64(0x94D049BB133111EB)
        x ^= x >> np.uint64(31)
    return x


def pair_hash(hr: np.ndarray, ht: np.ndarray) -> np.ndarray:
    """Seeded tie-break hash of (regulator, target) pairs from per-gene hashes."""
    with np.errstate(over="ignore"):
        return _mix64(hr ^ _mix64(ht + np.uint64(0x9E3779B97F4A7C15)))


def quantile_bins(values, n: int) -> np.ndarray:
    """1..n quantile bins from average ranks; tied values always share a bin."""
    v = np.asarray(values, dtype=float)
    if len(v) == 0:
        return np.zeros(0, dtype=np.int64)
    r = stats.rankdata(v, method="average")
    b = np.floor(n * (r - 0.5) / len(v)).astype(np.int64)
    return np.clip(b, 0, n - 1) + 1


def odds_ratio(a, b, c, d):
    """Sample OR ad/bc; if any cell is 0, Haldane-Anscombe (+0.5 on all cells). Vectorised."""
    a, b, c, d = (np.asarray(x, dtype=float) for x in (a, b, c, d))
    if (np.minimum(np.minimum(a, b), np.minimum(c, d)) < -1e-6).any():
        raise ValueError("negative 2x2 cell")
    a, b, c, d = (np.maximum(x, 0.0) for x in (a, b, c, d))
    zero = (a == 0) | (b == 0) | (c == 0) | (d == 0)
    h = np.where(zero, 0.5, 0.0)
    with np.errstate(all="ignore"):     # Accelerate leaves spurious FP flags (DESIGN §10)
        return ((a + h) * (d + h)) / ((b + h) * (c + h))


def log_or(a, b, c, d):
    o = odds_ratio(a, b, c, d)
    with np.errstate(all="ignore"):
        return np.log(o)


def _safe_log(x) -> float:
    with np.errstate(all="ignore"):
        return float(np.log(x))


def table_from_counts(a, s, p, n):
    """2x2 (a, b, c, d) from: a = selected & positive, s = selected, p = positive, n = all pairs."""
    return a, s - a, p - a, n - s - p + a


def fisher(a, b, c, d) -> dict:
    t = [[int(a), int(b)], [int(c), int(d)]]
    _, pg = stats.fisher_exact(t, alternative="greater")
    _, p2 = stats.fisher_exact(t, alternative="two-sided")
    return {"fisher_p_greater": float(pg), "fisher_p_two_sided": float(p2)}


def mantel_haenszel(a, b, c, d) -> dict:
    """MH pooled OR, Robins-Breslow-Greenland SE of log OR, CMH chi-square (no continuity)."""
    a, b, c, d = (np.asarray(x, dtype=float) for x in (a, b, c, d))
    n = a + b + c + d
    k = n > 0
    a, b, c, d, n = a[k], b[k], c[k], d[k], n[k]
    R = a * d / n
    S = b * c / n
    sR, sS = R.sum(), S.sum()
    out = {"cmh_or": float(sR / sS) if sS > 0 else float("inf"), "n_strata": int(k.sum())}
    if sR > 0 and sS > 0:
        P = (a + d) / n
        Q = (b + c) / n
        var = ((P * R).sum() / (2 * sR ** 2) + (P * S + Q * R).sum() / (2 * sR * sS)
               + (Q * S).sum() / (2 * sS ** 2))
        se = math.sqrt(var)
        out.update(cmh_log_se=se, cmh_ci95=[float(math.exp(math.log(sR / sS) - 1.959964 * se)),
                                            float(math.exp(math.log(sR / sS) + 1.959964 * se))])
    m = n > 1
    ea = ((a + b) * (a + c) / n)[m]
    va = ((a + b) * (c + d) * (a + c) * (b + d) / (n ** 2 * (n - 1)))[m]
    if va.sum() > 0:
        chi2 = (a[m].sum() - ea.sum()) ** 2 / va.sum()
        out.update(cmh_chi2=float(chi2), cmh_p=float(stats.chi2.sf(chi2, 1)))
    return out


def sha256_file(path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def git_hash(root: Path = ROOT) -> str:
    try:
        h = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True,
                           timeout=10).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=root,
                               capture_output=True, text=True, timeout=10).stdout.strip()
        return (h or "nogit") + ("+dirty" if dirty else "")
    except Exception:  # pragma: no cover
        return "nogit"


# =====================================================================================
# B. building the benchmark tables (pure functions)
# =====================================================================================

def pie_key_map(pie_map: pd.DataFrame) -> pd.Series:
    """PIE identifier -> DepMap gene key ('SYMBOL (ENTREZ)'); NaN if unmapped. Asserts 1:1."""
    key = pie_map["depmap_cn_column"].where(pie_map["depmap_cn_column"].notna(),
                                            pie_map["depmap_expr_column"])
    s = pd.Series(key.values, index=pie_map["pie_id"].values)
    if s.index.duplicated().any():
        raise ValueError("duplicated pie_id in pie_gene_map")
    m = s.dropna()
    if m.duplicated().any():
        raise ValueError(f"PIE ids map many-to-one to DepMap keys: {m[m.duplicated(keep=False)].head()}")
    return s


KEY_COLS = ("gene", "gene_key", "depmap_column", "key")


def normalise_flags(flags: pd.DataFrame) -> pd.DataFrame:
    """Index frozen flags by gene key; arm -> arm_id like '17q'; booleans coerced."""
    f = flags.copy()
    if f.index.name not in KEY_COLS:
        keycol = next((c for c in KEY_COLS if c in f.columns), None)
        if keycol is None:
            raise ValueError(f"frozen flags need a gene-key column, one of {KEY_COLS}")
        f = f.set_index(keycol)
    f.index.name = "gene"
    f["chrom"] = f["chrom"].astype("string").str.removeprefix("chr")
    arm = f["arm"].astype("string")
    short = arm.isin(["p", "q"])
    f["arm_id"] = arm.where(~short, f["chrom"] + arm)
    for c in ("variable_cn", "measured", "has_expr", "has_cn", "is_par", "is_y", "is_mt"):
        if c in f.columns:
            f[c] = f[c].map(lambda x: str(x).strip().lower() in ("true", "1", "yes")).astype(bool)
    return f


def _placed_ok(f: pd.DataFrame) -> pd.Series:
    """Autosome or X, not PAR, not Y, not MT, with an arm."""
    ok = f["chrom"].isin(CHROMS) & f["arm_id"].notna()
    for c in ("is_par", "is_y", "is_mt"):
        if c in f.columns:
            ok &= ~f[c].astype(bool)
    return ok.astype(bool)


def build_universe_regulators(pie_map: pd.DataFrame, flags: pd.DataFrame, line_targets: Iterable[str]
                              ) -> tuple[pd.DataFrame, dict]:
    """Regulators of one line's universe (DESIGN §6): PIE `target` in that line, 1:1 DepMap key,
    autosome/X, frozen variable-CN flag. Returns (table, exclusion counts)."""
    f = normalise_flags(flags)
    km = pie_key_map(pie_map)
    t = pd.Index(sorted(set(line_targets)))
    key = km.reindex(t)
    counts = {"pie_targets": len(t), "unmapped": int(key.isna().sum())}
    mapped = key.dropna()
    in_flags = mapped.isin(f.index)
    counts["mapped_not_in_flags"] = int((~in_flags).sum())
    mapped = mapped[in_flags]
    ff = f.loc[mapped.values]
    placed = _placed_ok(ff).values
    counts["excluded_Y_PAR_MT_unplaced"] = int((~placed).sum())
    var = ff["variable_cn"].values
    counts["not_variable_cn"] = int((placed & ~var).sum())
    keep = placed & var
    out = pd.DataFrame({"regulator": mapped.values[keep], "pie_id": mapped.index[keep],
                        "chrom": ff["chrom"].values[keep], "arm": ff["arm_id"].values[keep]})
    counts["universe_regulators"] = len(out)
    return out.sort_values("regulator").reset_index(drop=True), counts


def build_universe_targets(pie_map: pd.DataFrame, flags: pd.DataFrame, consistent_features: Iterable[str],
                           gene_pos: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Targets of one line's universe: consistently tested PIE feature (control mean >= 0.05), 1:1 DepMap
    key, frozen `measured` flag, autosome/X (no Y/PAR/MT). `gene_pos`: index gene key, columns start, end."""
    f = normalise_flags(flags)
    km = pie_key_map(pie_map)
    feats = pd.Index(sorted(set(consistent_features)))
    key = km.reindex(feats)
    counts = {"consistent_features": len(feats), "unmapped": int(key.isna().sum())}
    mapped = key.dropna()
    in_flags = mapped.isin(f.index)
    counts["mapped_not_in_flags"] = int((~in_flags).sum())
    mapped = mapped[in_flags]
    ff = f.loc[mapped.values]
    placed = _placed_ok(ff).values
    counts["excluded_Y_PAR_MT_unplaced"] = int((~placed).sum())
    meas = ff["measured"].values
    counts["not_measured"] = int((placed & ~meas).sum())
    keep = placed & meas
    out = pd.DataFrame({"target": mapped.values[keep], "pie_id": mapped.index[keep],
                        "chrom": ff["chrom"].values[keep], "arm": ff["arm_id"].values[keep]})
    pos = gene_pos.reindex(out["target"])
    if pos[["start", "end"]].isna().any().any():
        raise ValueError("universe targets without coordinates: "
                         f"{list(out['target'][pos['start'].isna().values][:5])}")
    out["mid"] = ((pos["start"].values + pos["end"].values) / 2).astype(np.int64)
    out["bin"] = out["chrom"] + ":" + (out["mid"] // BIN_BP).astype(str)
    counts["universe_targets"] = len(out)
    return out.sort_values("target").reset_index(drop=True), counts


def own_gene_class(own_tested, own_padj, own_lfc) -> np.ndarray:
    """(i) own feature tested and knocked down (p_adj < 0.05, LFC < 0); (ii) tested, not knocked down;
    (iii) own feature not tested (not expressed)."""
    tested = np.asarray(own_tested, dtype=bool)
    padj = np.asarray(own_padj, dtype=float)
    lfc = np.asarray(own_lfc, dtype=float)
    kd = tested & (padj < PADJ_POS) & (lfc < 0)
    return np.where(kd, "i", np.where(tested, "ii", "iii"))


def build_positives(de_sig: pd.DataFrame, regs: pd.DataFrame, tgts: pd.DataFrame) -> pd.DataFrame:
    """Sparse positives (DESIGN §6; INTERFACES 'Benchmark tables').

    de_sig: PIE rows with columns target, feature, log2_fold_change, p_adj, Ueffect; it must contain
    every row with p_adj < 0.05 (the cached `loose_*` table has all p_adj < 0.1 rows).
    Kept: trans pairs inside the universe with p_adj < 0.05 (every sensitivity set is a subset)."""
    rmap = dict(zip(regs["pie_id"], regs["regulator"]))
    tmap = dict(zip(tgts["pie_id"], tgts["target"]))
    d = de_sig[(de_sig["p_adj"] < PADJ_POS) & de_sig["target"].isin(rmap) & de_sig["feature"].isin(tmap)]
    out = pd.DataFrame({"regulator": d["target"].map(rmap).values, "target": d["feature"].map(tmap).values,
                        "p_adj": d["p_adj"].values.astype(float),
                        "lfc_sign": np.sign(d["log2_fold_change"].values).astype(np.int8),
                        "ueffect": d["Ueffect"].values.astype(float)})
    rc = dict(zip(regs["regulator"], regs["chrom"]))
    tc = dict(zip(tgts["target"], tgts["chrom"]))
    trans = out["regulator"].map(rc).values != out["target"].map(tc).values
    out = out[trans & (out["regulator"] != out["target"]).values].copy()
    out["pos_primary"] = out["p_adj"] < PADJ_POS
    out["pos_ueffect"] = out["pos_primary"] & (out["ueffect"].abs() >= UEFFECT_MIN)
    out["pos_p01"] = out["p_adj"] < PADJ_STRICT
    if out.duplicated(["regulator", "target"]).any():
        raise ValueError("duplicated (regulator, target) in positives")
    return out.sort_values(["regulator", "target"]).reset_index(drop=True)


def degree_deciles(regs: pd.DataFrame, positives: pd.DataFrame, split_col: str = "split") -> pd.Series:
    """Benchmark degree (# primary positives) deciles, computed within split (evaluation-only)."""
    deg = positives[positives["pos_primary"]].groupby("regulator").size()
    d = regs["regulator"].map(deg).fillna(0).astype(int)
    out = pd.Series(0, index=regs.index, dtype=np.int64)
    for _, idx in regs.groupby(split_col).groups.items():
        out.loc[idx] = quantile_bins(d.loc[idx].values, 10)
    return out


# =====================================================================================
# C. held-out split
# =====================================================================================

def make_split(arm_counts: pd.DataFrame, seed: int = SEED) -> pd.DataFrame:
    """Arm-level split (DESIGN §6, §8.2).

    arm_counts: columns arm, n_hct116, n_hek293t (universe regulators per arm).
    Units: arms of chr1-22/X with >= 1 universe regulator in either line, sorted by the number of HCT116
    universe regulators (descending; ties in genomic order). In consecutive triplets one arm is drawn
    uniformly to be TEST (np.random.default_rng(seed).integers(3) per triplet, in order). A final partial
    group of m < 3 arms uses the same draw; it gets a TEST arm only if the draw is < m.
    TRAIN arms alternate fold 1, 2, 1, 2, ... in the sorted order."""
    a = arm_counts.copy()
    a = a[(a["n_hct116"] + a["n_hek293t"]) > 0]
    a = a[a["arm"].map(lambda x: x[:-1] in CHROMS)]
    a["_g"] = a["arm"].map(arm_sort_key)
    a = a.sort_values(["n_hct116", "_g"], ascending=[False, True], kind="mergesort").drop(columns="_g")
    a = a.reset_index(drop=True)
    a["order"] = np.arange(len(a))
    a["triplet"] = a["order"] // 3
    rng = np.random.default_rng(seed)
    split = np.array(["TRAIN"] * len(a), dtype=object)
    draws = []
    for t in range(int(math.ceil(len(a) / 3))):
        k = int(rng.integers(3))
        draws.append(k)
        if 3 * t + k < len(a):
            split[3 * t + k] = "TEST"
    a["split"] = split
    fold = np.zeros(len(a), dtype=np.int64)
    tr = np.flatnonzero(split == "TRAIN")
    fold[tr] = np.arange(len(tr)) % 2 + 1
    a["fold"] = fold
    return a


# =====================================================================================
# D. loading benchmark tables and score files (with the TEST guard)
# =====================================================================================

class SealedSplitError(PermissionError):
    """Raised when TEST labels are requested outside scripts/final_eval.py."""


class _FinalEvalToken:
    __slots__ = ("label_dir", "mock")

    def __init__(self, label_dir: Path, mock: bool):
        self.label_dir = Path(label_dir)
        self.mock = mock


def final_eval_token(label_dir, mock: bool = False) -> _FinalEvalToken:
    """Only scripts/final_eval.py may obtain this token (checked on the running program's path)."""
    prog = Path(sys.argv[0]).resolve() if sys.argv and sys.argv[0] else Path("")
    if not (prog.name == "final_eval.py" and prog.parent.name == "scripts"):
        raise SealedSplitError("TEST evaluation is only allowed from scripts/final_eval.py")
    return _FinalEvalToken(label_dir, mock)


def bench_dir(root: Path = ROOT) -> Path:
    return Path(root) / "data" / "processed" / "bench"


def _check_token(split: str, token) -> None:
    if split != "TRAIN" and not isinstance(token, _FinalEvalToken):
        raise SealedSplitError(f"split={split!r} refused: TEST labels are sealed (DESIGN §6). "
                               "Only scripts/final_eval.py can evaluate TEST.")


def read_positives(path, token=None) -> pd.DataFrame:
    """Read a positives file. Any file whose name mentions TEST needs the final-eval token."""
    if "TEST" in Path(path).name.upper() and not isinstance(token, _FinalEvalToken):
        raise SealedSplitError(f"refusing to read {path}: TEST labels are sealed")
    return pd.read_parquet(path)


@dataclass
class Bench:
    line: str
    split: str
    regs: pd.DataFrame          # regulator, arm, chrom, own_class, cn_tercile, degree_decile, fold, ...
    tgts: pd.DataFrame          # target, chrom, bin, expr_decile, ...
    pos: pd.DataFrame           # regulator, target, p_adj, lfc_sign, ueffect, pos_* flags
    files: dict = field(default_factory=dict)   # path -> sha256 (for the log)


def load_bench(line: str, split: str = "TRAIN", root: Path = ROOT, token=None, mock: bool = False) -> Bench:
    """Load one line x split. TRAIN is open; TEST needs the final-eval token.
    mock=True (final_eval --mock): TRAIN fold-2 arms relabelled as a pseudo-TEST split."""
    if line not in LINES:
        raise ValueError(line)
    if mock:
        split = "TRAIN"
    _check_token(split, token)
    bd = bench_dir(root)
    rp = bd / f"universe_regulators_{line}.tsv"
    tp = bd / f"universe_targets_{line}.tsv"
    regs = pd.read_csv(rp, sep="\t", dtype={"chrom": str})
    tgts = pd.read_csv(tp, sep="\t", dtype={"chrom": str})
    regs = regs[regs["split"] == split].reset_index(drop=True)
    files = {str(rp.relative_to(root)): sha256_file(rp), str(tp.relative_to(root)): sha256_file(tp)}
    if split == "TRAIN":
        pp = bd / f"positives_{line}_TRAIN.parquet"
        pos = read_positives(pp)
        files[str(pp.relative_to(root))] = sha256_file(pp)
    else:
        ld = token.label_dir
        pp = ld / f"positives_{line}_TEST.parquet"
        lp = ld / f"regulator_labels_{line}_TEST.tsv"
        pos = read_positives(pp, token)
        lab = pd.read_csv(lp, sep="\t")
        regs = regs.drop(columns=[c for c in ("own_class", "degree_decile") if c in regs.columns])
        regs = regs.merge(lab, on="regulator", how="left", validate="1:1")
        if regs[["own_class", "degree_decile"]].isna().any().any():
            raise ValueError("TEST regulator labels incomplete")
        files[str(pp)] = sha256_file(pp)
        files[str(lp)] = sha256_file(lp)
    if mock:
        regs = regs[regs["fold"] == 2].reset_index(drop=True)
        regs["split"] = "PSEUDO_TEST"
        pos = pos[pos["regulator"].isin(regs["regulator"])].reset_index(drop=True)
        regs["degree_decile"] = degree_deciles(regs, pos).values
        split = "PSEUDO_TEST"
    return Bench(line=line, split=split, regs=regs, tgts=tgts, pos=pos, files=files)


def resolve_run_dir(run, root: Path = ROOT) -> Path:
    p = Path(run)
    if p.is_dir():
        return p
    q = Path(root) / "data" / "results" / "scores" / str(run)
    if q.is_dir():
        return q
    raise FileNotFoundError(f"score directory not found: {run}")


def load_scores(run_dir) -> dict:
    """Read scores.npz (+ meta.json) per INTERFACES 'Score files'."""
    run_dir = Path(run_dir)
    z = np.load(run_dir / "scores.npz", allow_pickle=False)
    out = {"regulators": np.asarray(z["regulators"]).astype(str),
           "targets": np.asarray(z["targets"]).astype(str),
           "score": z["score"],
           "s_cis": np.asarray(z["s_cis"]).astype(np.int8)}
    if "call" in z.files:
        out["call"] = z["call"].astype(bool)
    if not set(np.unique(out["s_cis"])) <= {-1, 1}:
        raise ValueError("s_cis must be in {-1, +1}")
    if out["score"].shape != (len(out["regulators"]), len(out["targets"])):
        raise ValueError("score shape does not match regulators x targets")
    mp = run_dir / "meta.json"
    out["meta"] = json.loads(mp.read_text()) if mp.exists() else {}
    out["run_dir"] = run_dir
    return out


def align_scores(sc: dict, regs: pd.Series, tgts: pd.Series, key: str = "score"):
    """Universe-aligned dense block of a [R, T] score-file array; absent genes -> NaN (False for masks)."""
    ri = pd.Index(sc["regulators"]).get_indexer(regs)
    ti = pd.Index(sc["targets"]).get_indexer(tgts)
    arr = sc[key]
    fill = False if arr.dtype == bool else np.nan
    out = np.full((len(ri), len(ti)), fill, dtype=arr.dtype if arr.dtype == bool else np.float32)
    rok, tok = ri >= 0, ti >= 0
    if rok.any() and tok.any():
        out[np.ix_(np.flatnonzero(rok), np.flatnonzero(tok))] = arr[np.ix_(ri[rok], ti[tok])]
    s_cis = np.ones(len(ri), dtype=np.int8)
    s_cis[rok] = sc["s_cis"][ri[rok]]
    return out, s_cis, int((~rok).sum()), int((~tok).sum())


# =====================================================================================
# E. evaluation engine
# =====================================================================================

class Evaluation:
    """All metrics for one score matrix on one (line, split) universe.

    Labels stay sparse (linear pair indices row * T + col). The only dense arrays are the score block and
    its ranking key (float32, one regulator split x universe targets)."""

    def __init__(self, regs: pd.DataFrame, tgts: pd.DataFrame, pos: pd.DataFrame, S: np.ndarray,
                 s_cis: np.ndarray, seed: int = SEED):
        self.regs = regs.reset_index(drop=True)
        self.tgts = tgts.reset_index(drop=True)
        self.seed = seed
        nr, nt = len(regs), len(tgts)
        if S.shape != (nr, nt):
            raise ValueError("score block shape mismatch")
        self.nr, self.nt = nr, nt
        chroms = pd.Index(sorted(set(self.regs["chrom"].astype(str)) | set(self.tgts["chrom"].astype(str)),
                                 key=chrom_order))
        self.chroms = chroms
        self.cr = chroms.get_indexer(self.regs["chrom"].astype(str))
        self.ct = chroms.get_indexer(self.tgts["chrom"].astype(str))
        self.nC = len(chroms)
        self.nT_c = np.bincount(self.ct, minlength=self.nC)
        self.n_i = nt - self.nT_c[self.cr]                      # universe pairs per regulator
        self.N = int(self.n_i.sum())
        self.S = S
        self.s_cis = np.asarray(s_cis, dtype=np.int8)
        key = np.abs(S).astype(np.float32)
        key[np.isnan(key)] = -1.0                               # NaN -> minimum
        for c in range(self.nC):
            r = np.flatnonzero(self.cr == c)
            t = np.flatnonzero(self.ct == c)
            if len(r) and len(t):
                key[np.ix_(r, t)] = -2.0                        # same chromosome (incl. A = B): masked
        self.key = key
        self.hr = key_hash(self.regs["regulator"], seed)
        self.ht = key_hash(self.tgts["target"], seed)
        # arms, bins, deciles
        self.arms = pd.Index(sorted(self.regs["arm"].unique(), key=arm_sort_key))
        self.ga = self.arms.get_indexer(self.regs["arm"])
        self.G = len(self.arms)
        self.bins = pd.Index(sorted(self.tgts["bin"].unique()))
        self.tb = self.bins.get_indexer(self.tgts["bin"])
        self.H = len(self.bins)
        self.bin_chrom = chroms.get_indexer(pd.Series(self.bins).str.split(":").str[0])
        self.nB = np.bincount(self.tb, minlength=self.H)
        self.te = (self.tgts["expr_decile"].astype(int).values - 1) if "expr_decile" in self.tgts else None
        # positives
        rix = pd.Index(self.regs["regulator"])
        tix = pd.Index(self.tgts["target"])
        r = rix.get_indexer(pos["regulator"])
        t = tix.get_indexer(pos["target"])
        ok = (r >= 0) & (t >= 0)
        p = pos[ok].copy()
        p["lin"] = r[ok].astype(np.int64) * nt + t[ok]
        p = p.sort_values("lin").reset_index(drop=True)
        if p["lin"].duplicated().any():
            raise ValueError("duplicate positives")
        if (self.key.ravel()[p["lin"].values] <= -2).any():
            raise ValueError("a positive pair is same-chromosome or A = B: positives must be trans")
        self.pos = p
        self._sel_cache: dict = {}

    # ---------- basic pieces ----------
    def n_universe_pairs(self) -> int:
        return self.N

    def topk(self, frac: float = PRIMARY_K) -> np.ndarray:
        """Sorted linear indices of the top-K pairs by |score| (NaN -> min; ties by seeded hash)."""
        if frac in self._sel_cache:
            return self._sel_cache[frac]
        K = max(1, int(round(frac * self.N)))
        flat = self.key.ravel()
        thr = np.partition(flat, flat.size - K)[flat.size - K]
        gt = np.flatnonzero(flat > thr)
        need = K - len(gt)
        eq = np.flatnonzero(flat == thr)
        if need < len(eq):
            h = pair_hash(self.hr[eq // self.nt], self.ht[eq % self.nt])
            eq = eq[np.argsort(h, kind="stable")[:need]]
        sel = np.sort(np.concatenate([gt, eq]))
        assert len(sel) == K and (flat[sel] > -2).all()
        self._sel_cache[frac] = sel
        return sel

    def pos_lin(self, flag: str = "pos_primary") -> np.ndarray:
        return self.pos.loc[self.pos[flag].astype(bool), "lin"].values

    def per_reg(self, sel: np.ndarray, pl: np.ndarray):
        """Per-regulator counts: selected, positive, both."""
        both = np.intersect1d(sel, pl, assume_unique=True)
        f = lambda x: np.bincount(x // self.nt, minlength=self.nr)
        return f(sel), f(pl), f(both), both

    def pooled(self, sel, pl, reg_mask=None) -> dict:
        s, p, a, _ = self.per_reg(sel, pl)
        m = slice(None) if reg_mask is None else reg_mask
        tab = table_from_counts(a[m].sum(), s[m].sum(), p[m].sum(), self.n_i[m].sum())
        return {"a": int(tab[0]), "b": int(tab[1]), "c": int(tab[2]), "d": int(tab[3]),
                "or": float(odds_ratio(*tab)), "log_or": float(log_or(*tab))}

    # ---------- permutation ----------
    def _strata_index(self, strata: np.ndarray):
        groups = pd.Series(np.arange(self.nr)).groupby(pd.Series(strata).astype(str).values).indices
        return [np.asarray(v) for v in groups.values()]

    def permutation(self, sel, pl, strata, n_perm: int = N_PERM, rng=None, perms=None,
                    return_null: bool = False) -> dict:
        """Regulator-label permutation (DESIGN §6). Label rows are reassigned among regulators within strata;
        each replicate uses the pairs trans for both A and pi(A); the observed OR of that replicate uses the
        same pair set; top-K selection fixed. p = (1 + #{OR_null >= OR_obs}) / (n_perm + 1)."""
        rng = rng if rng is not None else spawn_rng(SPAWN["perm_regulator"], self.seed)
        nr, nt, nC = self.nr, self.nt, self.nC
        sr, sc = sel // nt, sel % nt
        pr, pc = pl // nt, pl % nt
        s_rc = np.bincount(sr * nC + self.ct[sc], minlength=nr * nC).reshape(nr, nC)
        l_rc = np.bincount(pr * nC + self.ct[pc], minlength=nr * nC).reshape(nr, nC)
        both = np.intersect1d(sel, pl, assume_unique=True)
        o_rc = np.bincount((both // nt) * nC + self.ct[both % nt], minlength=nr * nC).reshape(nr, nC)
        s_tot, l_tot, o_tot = s_rc.sum(1), l_rc.sum(1), o_rc.sum(1)
        Sel = sp.csr_matrix((np.ones(len(sel), np.int64), (sr, sc)), shape=(nr, nt))
        Pos = sp.csr_matrix((np.ones(len(pl), np.int64), (pr, pc)), shape=(nr, nt))
        M = (Sel @ Pos.T).toarray()                     # M[i, j] = |sel(i) ∩ pos(j)|
        groups = self._strata_index(strata)
        i = np.arange(nr)
        ci = self.cr
        full = table_from_counts(o_tot.sum(), s_tot.sum(), l_tot.sum(), self.N)
        or_full = float(odds_ratio(*full))
        if perms is not None:                           # explicit permutations (tests)
            n_perm = len(perms)
        null = np.empty(n_perm)
        obs = np.empty(n_perm)
        tabs = []
        for k in range(n_perm):
            if perms is not None:
                j = np.asarray(perms[k], dtype=np.int64)
            else:
                j = np.empty(nr, dtype=np.int64)
                for g in groups:
                    j[g] = g[rng.permutation(len(g))]
            cj = self.cr[j]
            n = (nt - self.nT_c[ci] - np.where(cj == ci, 0, self.nT_c[cj])).sum()
            s = (s_tot - s_rc[i, cj]).sum()
            pn = (l_tot[j] - l_rc[j, ci]).sum()
            an = M[i, j].sum()
            po = (l_tot - l_rc[i, cj]).sum()
            ao = (o_tot - o_rc[i, cj]).sum()
            null[k] = odds_ratio(*table_from_counts(an, s, pn, n))
            obs[k] = odds_ratio(*table_from_counts(ao, s, po, n))
            if return_null:
                tabs.append((table_from_counts(an, s, pn, n), table_from_counts(ao, s, po, n)))
        ge = int(np.sum(null >= obs * (1 - 1e-12)))
        med = float(np.median(null))
        out = {"perm_p": (1 + ge) / (n_perm + 1), "n_perm": n_perm, "or": or_full,
               "median_null_or": med, "or_over_median_null": or_full / med if med > 0 else float("inf"),
               "null_or_q05_q95": [float(np.quantile(null, 0.05)), float(np.quantile(null, 0.95))]}
        if return_null:
            out.update(null=null, obs=obs, tables=tabs)
        return out

    # ---------- jackknife / folds / bootstrap ----------
    def arm_counts(self, sel, pl):
        s, p, a, _ = self.per_reg(sel, pl)
        f = lambda x: np.bincount(self.ga, weights=x, minlength=self.G)
        return f(a), f(s), f(p), f(self.n_i.astype(float))

    def jackknife(self, sels: list, pls: list) -> dict:
        """Delete-one-arm jackknife of log OR (one or two score files -> paired difference) with t_{G-1}."""
        G = self.G
        th_full, th_del = [], []
        for sel, pl in zip(sels, pls):
            a, s, p, n = self.arm_counts(sel, pl)
            th_full.append(log_or(*table_from_counts(a.sum(), s.sum(), p.sum(), n.sum())))
            th_del.append(log_or(*table_from_counts(a.sum() - a, s.sum() - s, p.sum() - p, n.sum() - n)))
        th = th_full[0] if len(sels) == 1 else th_full[0] - th_full[1]
        td = th_del[0] if len(sels) == 1 else th_del[0] - th_del[1]
        if G < 2:
            return {"G": G, "estimate": float(th)}
        se = math.sqrt((G - 1) / G * np.sum((td - td.mean()) ** 2))
        tq = stats.t.ppf(0.975, G - 1)
        tq90 = stats.t.ppf(0.95, G - 1)
        return {"G": G, "estimate": float(th), "se": se, "ci95": [float(th - tq * se), float(th + tq * se)],
                "lb90_one_sided": float(th - tq90 * se), "t_crit": float(tq)}

    def fold_log_or(self, sel, pl) -> dict:
        if "fold" not in self.regs or self.regs["fold"].isna().all():
            return {}
        out = {}
        s, p, a, _ = self.per_reg(sel, pl)
        for f in (1, 2):
            m = (self.regs["fold"] == f).values
            out[f"fold{f}"] = float(log_or(*table_from_counts(a[m].sum(), s[m].sum(), p[m].sum(),
                                                              self.n_i[m].sum())))
        return out

    def arm_bin_counts(self, sel, pl):
        """(G x H) counts of selected, positive, both and all universe pairs."""
        G, H, nt = self.G, self.H, self.nt
        f = lambda lin: np.bincount(self.ga[lin // nt] * H + self.tb[lin % nt], minlength=G * H
                                    ).reshape(G, H).astype(float)
        both = np.intersect1d(sel, pl, assume_unique=True)
        regs_per_arm = np.bincount(self.ga, minlength=G)
        arm_chrom = np.zeros(G, dtype=np.int64)
        arm_chrom[self.ga] = self.cr
        n = regs_per_arm[:, None] * self.nB[None, :] * (arm_chrom[:, None] != self.bin_chrom[None, :])
        return f(both), f(sel), f(pl), n.astype(float)

    def boot_weights(self, n_boot: int = N_BOOT, rng=None):
        rng = rng if rng is not None else spawn_rng(SPAWN["bootstrap"], self.seed)
        Wa = rng.multinomial(self.G, np.full(self.G, 1 / self.G), size=n_boot).astype(float)
        Wb = rng.multinomial(self.H, np.full(self.H, 1 / self.H), size=n_boot).astype(float)
        return Wa, Wb

    def bootstrap_log_or(self, sel, pl, W) -> np.ndarray:
        """Two-way bootstrap (regulator arms x target 5 Mb bins) replicates of log OR."""
        Wa, Wb = W
        with np.errstate(all="ignore"):     # spurious Accelerate matmul flags (DESIGN §10)
            q = [((Wa @ Q) * Wb).sum(1) for Q in self.arm_bin_counts(sel, pl)]
        if not all(np.isfinite(x).all() for x in q):
            raise FloatingPointError("non-finite bootstrap counts")
        return log_or(*table_from_counts(*q))

    # ---------- CMH ----------
    def cmh_tables(self, sel, pl, by: str = "regulator"):
        """Per-stratum 2x2: strata = regulator x target-expression decile ('regulator') or
        own-gene class x target-expression decile ('own_class')."""
        nt, nr, E = self.nt, self.nr, 10
        f = lambda lin: np.bincount((lin // nt) * E + self.te[lin % nt], minlength=nr * E).reshape(nr, E)
        both = np.intersect1d(sel, pl, assume_unique=True)
        a, s, p = f(both), f(sel), f(pl)
        nE_c = np.zeros((E, self.nC), dtype=np.int64)
        np.add.at(nE_c, (self.te, self.ct), 1)
        n = nE_c.sum(1)[None, :] - nE_c[:, self.cr].T          # (nr, E)
        if by == "own_class":
            cls = pd.Index(sorted(self.regs["own_class"].astype(str).unique()))
            gi = cls.get_indexer(self.regs["own_class"].astype(str))
            agg = lambda x: np.vstack([x[gi == g].sum(0) for g in range(len(cls))])
            a, s, p, n = agg(a), agg(s), agg(p), agg(n)
        return table_from_counts(a.ravel(), s.ravel(), p.ravel(), n.ravel())

    def cmh(self, sel, pl, by: str = "regulator") -> dict:
        return mantel_haenszel(*self.cmh_tables(sel, pl, by))

    def bootstrap_log_cmh(self, sel, pl, W, chunk: int = 100) -> np.ndarray:
        """Two-way bootstrap replicates of the regulator x expression-decile CMH log OR."""
        Wa, Wb = W
        B = Wa.shape[0]
        nt, E, H = self.nt, 10, self.H
        both = np.intersect1d(sel, pl, assume_unique=True)
        sk = lambda lin: (lin // nt) * E + self.te[lin % nt]
        active = np.unique(np.concatenate([sk(sel), sk(pl)]))
        kix = pd.Index(active)
        ki = active // E
        ke = active % E
        mk = lambda lin: sp.csr_matrix((np.ones(len(lin)), (kix.get_indexer(sk(lin)), self.tb[lin % nt])),
                                       shape=(len(active), H))
        Qa, Qs, Qp = mk(both), mk(sel), mk(pl)
        nEB = np.zeros((E, H))
        np.add.at(nEB, (self.te, self.tb), 1)
        with np.errstate(all="ignore"):     # spurious Accelerate matmul flags (DESIGN §10)
            NE = Wb @ nEB.T                                     # (B, E)
            NEC = np.zeros((B, E, self.nC))
            for c in range(self.nC):
                hc = self.bin_chrom == c
                NEC[:, :, c] = Wb[:, hc] @ nEB[:, hc].T
        if not (np.isfinite(NE).all() and np.isfinite(NEC).all()):
            raise FloatingPointError("non-finite bootstrap counts")
        out = np.empty(B)
        for b0 in range(0, B, chunk):
            b1 = min(B, b0 + chunk)
            wb = Wb[b0:b1].T                                    # (H, b)
            with np.errstate(all="ignore"):
                a = (Qa @ wb).T                                 # (b, K)
                s = (Qs @ wb).T
                p = (Qp @ wb).T
            n = NE[b0:b1][:, ke] - NEC[b0:b1][:, ke, self.cr[ki]]
            A, Bc, C, D = table_from_counts(a, s, p, n)
            with np.errstate(divide="ignore", invalid="ignore"):
                num = np.where(n > 0, A * D / n, 0.0)
                den = np.where(n > 0, Bc * C / n, 0.0)
            w = Wa[b0:b1][:, self.ga[ki]]                       # (b, K) arm multiplicities
            with np.errstate(all="ignore"):
                out[b0:b1] = np.log((w * num).sum(1) / (w * den).sum(1))
        return out

    # ---------- other secondary metrics ----------
    def auroc_per_regulator(self, pl, min_pos: int = 10) -> dict:
        p = np.bincount(pl // self.nt, minlength=self.nr)
        rows = np.flatnonzero(p >= min_pos)
        aucs = []
        for i in rows:
            valid = self.ct != self.cr[i]
            k = self.key[i, valid]
            y = np.zeros(self.nt, dtype=bool)
            y[pl[pl // self.nt == i] % self.nt] = True
            y = y[valid]
            r = stats.rankdata(k)
            n1, n0 = y.sum(), (~y).sum()
            aucs.append((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
        return {"mean_auroc": float(np.mean(aucs)) if aucs else float("nan"), "n_regulators": len(rows)}

    def top10_share(self, sel) -> float:
        s = np.bincount(sel // self.nt, minlength=self.nr)
        return float(np.sort(s)[::-1][:10].sum() / len(sel))

    def pred_kd_sign(self, lin) -> np.ndarray:
        """Predicted knockdown sign = -sign(score) * s_cis (0 where the score is NaN or 0)."""
        v = self.S.ravel()[lin]
        sg = np.where(np.isnan(v), 0, np.sign(v)).astype(np.int8)
        return (-sg * self.s_cis[lin // self.nt]).astype(np.int8)

    def sign_metrics(self, sel, flag: str = "pos_primary", n_perm: int = N_PERM, rng=None,
                     min_overlap: int = 5) -> dict:
        rng = rng if rng is not None else spawn_rng(SPAWN["perm_sign"], self.seed)
        pp = self.pos[self.pos[flag].astype(bool)]
        pl = pp["lin"].values
        lfc = pp["lfc_sign"].values.astype(np.int8)
        pred_pos = self.pred_kd_sign(pl)
        conc_pl = pl[pred_pos == lfc]
        s, p, a, _ = self.per_reg(sel, conc_pl)
        tab = table_from_counts(a.sum(), s.sum(), p.sum(), self.N)
        out = {"sign_concordant_or": float(odds_ratio(*tab)), "sign_concordant_table": [int(x) for x in tab],
               **{f"sign_concordant_{k}": v for k, v in fisher(*tab).items()}}
        ov = np.isin(pl, sel, assume_unique=True)
        pred, lf = pred_pos[ov], lfc[ov]
        reg = pl[ov] // self.nt
        n = len(pred)
        out["n_overlap"] = int(n)
        if n == 0:
            return out
        conc = pred == lf
        q = float(np.mean(pred == -1))
        u = float(np.mean(lf == -1))
        c0 = q * u + (1 - q) * (1 - u)
        dd = int(np.sum((pred == -1) & (lf == -1)))
        du = int(np.sum((pred == -1) & (lf == 1)))
        ud = int(np.sum((pred == 1) & (lf == -1)))
        uu = int(np.sum((pred == 1) & (lf == 1)))
        obs = int(conc.sum())
        null = np.empty(n_perm)
        for k in range(n_perm):
            order = np.lexsort((rng.random(n), reg))            # reg is sorted -> within-regulator shuffle
            null[k] = np.sum(pred == lf[order])
        rr = pd.Series(conc).groupby(reg)
        cnt = rr.size()
        rw = rr.mean()[cnt >= min_overlap]
        out.update(concordance=float(conc.mean()), q_pred_down=q, u_lfc_down=u, c0=c0,
                   concordance_minus_c0=float(conc.mean() - c0),
                   sign_2x2=[dd, du, ud, uu], sign_or=float(odds_ratio(dd, du, ud, uu)),
                   sign_perm_p=(1 + int(np.sum(null >= obs))) / (n_perm + 1), sign_perm_null_mean=float(null.mean() / n),
                   regulator_weighted_concordance=float(rw.mean()) if len(rw) else float("nan"),
                   n_regulators_weighted=int(len(rw)))
        return out

    def callset(self, call_block: np.ndarray, flag: str = "pos_primary") -> dict:
        """OR of a boolean call mask (universe-aligned [nr, nt]); same-chromosome entries ignored."""
        valid = self.key > -2
        lin = np.flatnonzero((call_block & valid).ravel())
        s, p, a, _ = self.per_reg(lin, self.pos_lin(flag))
        tab = table_from_counts(a.sum(), s.sum(), p.sum(), self.N)
        return {"n_calls": int(len(lin)), "or": float(odds_ratio(*tab)), "table": [int(x) for x in tab],
                **fisher(*tab)}

    def subset(self, reg_mask=None, tgt_mask=None) -> "Evaluation":
        rm = np.ones(self.nr, bool) if reg_mask is None else np.asarray(reg_mask, bool)
        tm = np.ones(self.nt, bool) if tgt_mask is None else np.asarray(tgt_mask, bool)
        regs = self.regs[rm]
        tgts = self.tgts[tm]
        pos = self.pos[self.pos["regulator"].isin(regs["regulator"]) & self.pos["target"].isin(tgts["target"])]
        return Evaluation(regs, tgts, pos.drop(columns="lin"), self.S[np.ix_(rm, tm)], self.s_cis[rm], self.seed)

    def hub_excluded(self, frac_hub: float = 0.01, k: float = PRIMARY_K) -> dict:
        """Drop the top 1 % of regulators and of targets by benchmark degree (evaluated split), re-rank."""
        pl = self.pos_lin("pos_primary")
        dr = np.bincount(pl // self.nt, minlength=self.nr)
        dt = np.bincount(pl % self.nt, minlength=self.nt)
        nr_drop = int(math.ceil(frac_hub * self.nr))
        nt_drop = int(math.ceil(frac_hub * self.nt))
        ro = np.lexsort((self.hr, -dr))[:nr_drop]
        to = np.lexsort((self.ht, -dt))[:nt_drop]
        rm = np.ones(self.nr, bool)
        rm[ro] = False
        tm = np.ones(self.nt, bool)
        tm[to] = False
        sub = self.subset(rm, tm)
        sel = sub.topk(k)
        res = sub.pooled(sel, sub.pos_lin("pos_primary"))
        res.update(fisher(res["a"], res["b"], res["c"], res["d"]), n_reg_dropped=nr_drop, n_tgt_dropped=nt_drop)
        return res


# =====================================================================================
# F. evaluate() / compare() with the tuning log
# =====================================================================================

LOG_COLUMNS = ["timestamp", "run_id", "variant", "kind", "line", "split", "params_hash", "git_hash",
               "data_checksums", "seed", "spawn", "n_perm", "n_boot", "metrics", "wall_s"]


def _append_log(row: dict, root: Path) -> None:
    lp = Path(root) / "results" / "tuning_log.tsv"
    lp.parent.mkdir(parents=True, exist_ok=True)
    new = not lp.exists()
    with open(lp, "a") as fh:
        if new:
            fh.write("\t".join(LOG_COLUMNS) + "\n")
        fh.write("\t".join(str(row.get(c, "")).replace("\t", " ").replace("\n", " ") for c in LOG_COLUMNS) + "\n")


def _params_hash(meta: dict, extra: dict) -> str:
    blob = json.dumps({"params": meta.get("params", {}), "eval": extra}, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    return x


def build_evaluation(sc: dict, bench: Bench, reg_subset=None, seed: int = SEED):
    regs = bench.regs if reg_subset is None else bench.regs[np.asarray(reg_subset, bool)]
    S, s_cis, miss_r, miss_t = align_scores(sc, regs["regulator"], bench.tgts["target"])
    pos = bench.pos[bench.pos["regulator"].isin(regs["regulator"])]
    ev = Evaluation(regs, bench.tgts, pos, S, s_cis, seed)
    return ev, {"n_regulators_missing_in_scores": miss_r, "n_targets_missing_in_scores": miss_t}


def compute_metrics(ev: Evaluation, *, n_perm: int = N_PERM, n_boot: int = N_BOOT, call_block=None,
                    full: bool = True) -> dict:
    m: dict = {"n_regulators": ev.nr, "n_targets": ev.nt, "n_pairs": ev.N, "n_arms": ev.G,
               "n_target_bins": ev.H}
    pl = ev.pos_lin("pos_primary")
    m["n_positives"] = {f: int(ev.pos[f].sum()) for f in POS_FLAGS}
    for k in K_FRACS:
        sel = ev.topk(k)
        r = ev.pooled(sel, pl)
        r.update(fisher(r["a"], r["b"], r["c"], r["d"]), K=len(sel))
        m[f"top{k:g}"] = r
    sel = ev.topk(PRIMARY_K)
    m["primary_or"] = m[f"top{PRIMARY_K:g}"]["or"]
    m["sensitivity_positives"] = {f: ev.pooled(sel, ev.pos_lin(f)) for f in ("pos_ueffect", "pos_p01")}
    m["folds_log_or"] = ev.fold_log_or(sel, pl)
    m["jackknife_log_or"] = ev.jackknife([sel], [pl])
    if call_block is not None:
        m["callset"] = ev.callset(call_block)
    if not full:
        return m
    strata = ev.regs["own_class"].astype(str) + "|" + ev.regs["cn_tercile"].astype(str)
    m["perm_regulator"] = ev.permutation(sel, pl, strata.values, n_perm,
                                         spawn_rng(SPAWN["perm_regulator"], ev.seed))
    strata2 = ev.regs["own_class"].astype(str) + "|" + ev.regs["degree_decile"].astype(str)
    m["perm_degree"] = ev.permutation(sel, pl, strata2.values, n_perm, spawn_rng(SPAWN["perm_degree"], ev.seed))
    W = ev.boot_weights(n_boot)
    bl = ev.bootstrap_log_or(sel, pl, W)
    m["bootstrap_log_or"] = {"ci95": [float(np.quantile(bl, .025)), float(np.quantile(bl, .975))],
                             "n_boot": n_boot}
    m["cmh_regulator_x_exprdecile"] = ev.cmh(sel, pl, "regulator")
    m["cmh_ownclass_x_exprdecile"] = ev.cmh(sel, pl, "own_class")
    m["auroc"] = ev.auroc_per_regulator(pl)
    m["top10_regulator_share"] = ev.top10_share(sel)
    m["hub_excluded"] = ev.hub_excluded()
    m["sign"] = ev.sign_metrics(sel, "pos_primary", n_perm)
    for cls in ("i", "iii"):
        msk = (ev.regs["own_class"].astype(str) == cls).values
        if msk.sum() >= 1:
            sub = ev.subset(msk)
            ss = sub.topk(PRIMARY_K)
            r = sub.pooled(ss, sub.pos_lin("pos_primary"))
            r.update(fisher(r["a"], r["b"], r["c"], r["d"]), n_regulators=int(msk.sum()))
            m[f"class_{cls}"] = r
    return m


def evaluate(run, line: str = PRIMARY_LINE, split: str = "TRAIN", *, variant: Optional[str] = None,
             n_perm: int = N_PERM, n_boot: int = N_BOOT, full: bool = True, reg_subset=None,
             subset_name: Optional[str] = None, call_mask: Optional[np.ndarray] = None, root: Path = ROOT,
             log: bool = True, write_json: bool = True, _token=None, _mock: bool = False) -> dict:
    """Evaluate one score file on one line x split; append to results/tuning_log.tsv.

    TEST is refused unless `_token` comes from scripts/final_eval.py (final_eval_token)."""
    t0 = time.time()
    root = Path(root)
    _check_token("TRAIN" if _mock else split, _token)
    run_dir = resolve_run_dir(run, root)
    bench = load_bench(line, split, root, token=_token, mock=_mock)
    sc = load_scores(run_dir)
    if call_mask is not None:
        sc["call"] = np.asarray(call_mask, bool)
    ev, miss = build_evaluation(sc, bench, reg_subset)
    call_block = align_scores(sc, ev.regs["regulator"], ev.tgts["target"], "call")[0] if "call" in sc else None
    m = compute_metrics(ev, n_perm=n_perm, n_boot=n_boot, call_block=call_block, full=full)
    m.update(miss)
    wall = time.time() - t0
    meta = sc["meta"]
    run_id = meta.get("run_id", run_dir.name)
    variant = variant or meta.get("variant", meta.get("variant_id", ""))
    evp = {"k_fracs": K_FRACS, "n_perm": n_perm, "n_boot": n_boot, "full": full, "subset": subset_name}
    checks = dict(bench.files)
    checks[str(run_dir / "scores.npz")] = sha256_file(run_dir / "scores.npz")
    out = {"run_id": run_id, "variant": variant, "line": line, "split": bench.split, "metrics": _jsonable(m),
           "wall_s": wall, "score_wall_s": meta.get("wall_time"), "threads": meta.get("threads"),
           "seed": SEED, "spawn": SPAWN, "params_hash": _params_hash(meta, evp), "git_hash": git_hash(root),
           "data_checksums": checks, "subset": subset_name}
    if write_json:
        od = root / "results" / ("final" if bench.split == "TEST" else "eval") / str(run_id)
        if _mock:
            od = root / "results" / "final_mock" / str(run_id)
        od.mkdir(parents=True, exist_ok=True)
        name = f"{line}_{bench.split}" + (f"_{subset_name}" if subset_name else "") + ".json"
        (od / name).write_text(json.dumps(_jsonable(out), indent=1))
    if log:
        key = ["n_pairs", "primary_or"]
        summ = {k: m[k] for k in key}
        summ.update({f"top{k:g}_or": m[f"top{k:g}"]["or"] for k in K_FRACS})
        summ["fisher_p"] = m[f"top{PRIMARY_K:g}"]["fisher_p_greater"]
        summ.update(m.get("folds_log_or", {}))
        if full:
            summ.update(perm_p=m["perm_regulator"]["perm_p"], perm_degree_p=m["perm_degree"]["perm_p"],
                        cmh_or=m["cmh_regulator_x_exprdecile"]["cmh_or"],
                        sign_concordant_or=m["sign"]["sign_concordant_or"],
                        concordance=m["sign"].get("concordance"), c0=m["sign"].get("c0"),
                        jk_ci95=m["jackknife_log_or"].get("ci95"))
        _append_log({"timestamp": _dt.datetime.now().isoformat(timespec="seconds"), "run_id": run_id,
                     "variant": variant, "kind": "evaluate" + (f":{subset_name}" if subset_name else ""),
                     "line": line, "split": bench.split, "params_hash": out["params_hash"],
                     "git_hash": out["git_hash"],
                     "data_checksums": json.dumps({Path(k).name: v[:16] for k, v in checks.items()}),
                     "seed": SEED, "spawn": json.dumps(SPAWN), "n_perm": n_perm, "n_boot": n_boot,
                     "metrics": json.dumps(_jsonable(summ)), "wall_s": f"{wall:.2f}"}, root)
    return out


def compare(run_candidate, run_incumbent, line: str = PRIMARY_LINE, split: str = "TRAIN", *,
            n_boot: int = N_BOOT, root: Path = ROOT, log: bool = True, _token=None, _mock: bool = False) -> dict:
    """Paired comparison of two score files on the same universe (DESIGN §6 CIs, §8.2 adoption rule):
    Δlog OR per TRAIN fold, paired arm jackknife (t_{G-1}), paired two-way bootstrap of Δlog OR (90 % lower
    bound) and of Δlog CMH OR (regulator x expression decile; 90 % upper bound)."""
    t0 = time.time()
    root = Path(root)
    _check_token("TRAIN" if _mock else split, _token)
    bench = load_bench(line, split, root, token=_token, mock=_mock)
    split = bench.split
    evs, ids = [], []
    for run in (run_candidate, run_incumbent):
        rd = resolve_run_dir(run, root)
        sc = load_scores(rd)
        ev, _ = build_evaluation(sc, bench)
        evs.append(ev)
        ids.append(sc["meta"].get("run_id", rd.name))
    sels = [e.topk(PRIMARY_K) for e in evs]
    pl = evs[0].pos_lin("pos_primary")
    lo = [e.pooled(s, pl)["log_or"] for e, s in zip(evs, sels)]
    folds = [e.fold_log_or(s, pl) for e, s in zip(evs, sels)]
    dfold = {k: folds[0][k] - folds[1][k] for k in folds[0]}
    jk = evs[0].jackknife(sels, [pl, pl])
    W = evs[0].boot_weights(n_boot)
    bd = evs[0].bootstrap_log_or(sels[0], pl, W) - evs[1].bootstrap_log_or(sels[1], pl, W)
    cm = [e.cmh(s, pl)["cmh_or"] for e, s in zip(evs, sels)]
    with np.errstate(invalid="ignore"):
        bc = evs[0].bootstrap_log_cmh(sels[0], pl, W) - evs[1].bootstrap_log_cmh(sels[1], pl, W)
    res = {"candidate": ids[0], "incumbent": ids[1], "line": line, "split": split,
           "delta_log_or": lo[0] - lo[1], "delta_log_or_folds": dfold, "paired_jackknife": jk,
           "boot_delta_log_or_lb90": float(np.quantile(bd, 0.05)),
           "boot_delta_log_or_ci95": [float(np.quantile(bd, .025)), float(np.quantile(bd, .975))],
           "delta_log_cmh": _safe_log(cm[0]) - _safe_log(cm[1]),
           "boot_delta_log_cmh_ub90": float(np.quantile(bc, 0.95)), "n_boot": n_boot}
    res["adopt_rule_8_2"] = bool(len(dfold) == 2 and all(v > 0 for v in dfold.values())
                                 and res["boot_delta_log_or_lb90"] > 0 and res["boot_delta_log_cmh_ub90"] > 0)
    res["wall_s"] = time.time() - t0
    if log:
        _append_log({"timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
                     "run_id": f"{ids[0]}|{ids[1]}", "variant": "", "kind": "compare", "line": line,
                     "split": split, "params_hash": "", "git_hash": git_hash(root),
                     "data_checksums": json.dumps({Path(k).name: v[:16] for k, v in bench.files.items()}),
                     "seed": SEED, "spawn": json.dumps(SPAWN), "n_perm": 0, "n_boot": n_boot,
                     "metrics": json.dumps(_jsonable({k: v for k, v in res.items() if k != "paired_jackknife"})),
                     "wall_s": f"{res['wall_s']:.2f}"}, root)
    return _jsonable(res)
