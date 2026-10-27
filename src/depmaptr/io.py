"""Data loading for the model engine (DESIGN §2, §3, §8.6).

A `Dataset` holds everything the engine needs, aligned to panel row order:
- `panel`: one row per sample (raw covariates, panel_s2.tsv schema);
- `relcn`: linear relative CN (n x G_cn), `cn_keys`;
- `expr`: log2(TPM+1) (n x G_e), `expr_keys`;
- `genes`: annotation indexed by gene key (chrom, arm, start, flags);
- `context_expr`: {context line name: Series of expression over expr_keys} for G2 / B4e.

The same structure is produced by `simulate.py`, so every downstream function runs on SIM and on
real data unchanged.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data" / "processed"
STARTING = ROOT / "starting_data"

CHROMS = [str(i) for i in range(1, 23)] + ["X"]
CHROM_RANK = {c: i for i, c in enumerate(CHROMS)}

# DESIGN §3 default CN transform: x = log2(clip(relCN, 2^-2, 2^3)).
CN_FLOOR = 2.0 ** -2
CN_CAP = 2.0 ** 3
HOMDEL_RELCN = 0.25          # DESIGN §4: deletion indicator 1(relCN_B < 0.25)
HOMDEL_MIN_LINES = 3         # indicator only when >= 3 lines carry it
ARM_MIN_GENES = 20           # DESIGN §4: arms with >= 20 genes

# DESIGN §5.3: context lines for G2 / B4e, asserted against Model.csv CellLineName.
CONTEXT_LINES = {"HCT116": ("ACH-000971", "HCT 116"), "HEK293T": ("ACH-001085", "HEK 293")}


@dataclass
class Dataset:
    panel: pd.DataFrame
    relcn: np.ndarray
    cn_keys: np.ndarray
    expr: np.ndarray
    expr_keys: np.ndarray
    genes: pd.DataFrame
    regulators: np.ndarray
    targets: np.ndarray
    context_expr: dict = field(default_factory=dict)
    sample_set: str = "S2"
    input_checksums: dict = field(default_factory=dict)
    sample_notes: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return self.relcn.shape[0]

    @property
    def model_ids(self) -> np.ndarray:
        return self.panel["ModelID"].to_numpy()

    def subset_rows(self, mask: np.ndarray, name: str) -> "Dataset":
        mask = np.asarray(mask, bool)
        return Dataset(panel=self.panel.loc[mask].reset_index(drop=True), relcn=self.relcn[mask],
                       cn_keys=self.cn_keys, expr=self.expr[mask], expr_keys=self.expr_keys, genes=self.genes,
                       regulators=self.regulators, targets=self.targets, context_expr=self.context_expr,
                       sample_set=name, input_checksums=self.input_checksums, sample_notes=dict(self.sample_notes))


# ----------------------------------------------------------------------------- transforms
def cn_transform(relcn: np.ndarray, floor: float = CN_FLOOR, cap: float = CN_CAP) -> np.ndarray:
    """x = log2(clip(relCN, floor, cap)) (DESIGN §3), float64."""
    return np.log2(np.clip(np.asarray(relcn, np.float64), floor, cap))


def rint(M: np.ndarray) -> np.ndarray:
    """Rank-based inverse-normal transform per column: Phi^-1((rank - 0.5)/n), average ranks for ties."""
    from scipy.stats import rankdata
    M = np.asarray(M, np.float64)
    n = M.shape[0]
    R = rankdata(M, axis=0, method="average")
    return special.ndtri((R - 0.5) / n)


def gene_order(genes: pd.DataFrame, keys) -> np.ndarray:
    """Sort gene keys by (chromosome, arm p<q, start); genes without a start go last within their arm."""
    g = genes.loc[list(keys)]
    o = pd.DataFrame({"k": list(keys), "c": g["chrom"].map(CHROM_RANK).to_numpy(),
                      "a": (g["arm"] == "q").astype(int).to_numpy(),
                      "s": g["start"].fillna(np.inf).to_numpy()})
    return o.sort_values(["c", "a", "s", "k"])["k"].to_numpy()


# ----------------------------------------------------------------------------- real data
def sha256(path: Path, block: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(block):
            h.update(chunk)
    return h.hexdigest()


def context_model_check(model_csv: Path = STARTING / "Model.csv") -> None:
    """DESIGN §5.3: the context ACH IDs must be HCT 116 / HEK 293 in Model.csv."""
    m = pd.read_csv(model_csv, usecols=["ModelID", "CellLineName", "StrippedCellLineName"], low_memory=False)
    m = m.set_index("ModelID")
    for name, (ach, cell_line) in CONTEXT_LINES.items():
        got = m.loc[ach, "CellLineName"]
        # Model.csv spells HEK293 as "HEK-293"; compare on alphanumerics only.
        norm = lambda s: "".join(ch for ch in str(s).upper() if ch.isalnum())
        if norm(got) != norm(cell_line):
            raise AssertionError(f"context {name}: {ach} is {got!r} in Model.csv, expected {cell_line!r}")


def load_real(sample_set: str = "S2", proc: Path = PROC, checksums: bool = True) -> Dataset:
    """Load panel_s2 + matched CN/expression, aligned to panel order; S1 / leave-HCT116-out are subsets."""
    files = {k: proc / v for k, v in {
        "panel": "panel_s2.tsv", "flags": "frozen_flags.tsv", "annotation": "gene_annotation.tsv",
        "cn": "depmap_cn_matched.parquet", "expr": "depmap_expr_matched.parquet"}.items()}
    panel = pd.read_csv(files["panel"], sep="\t", low_memory=False)
    flags = pd.read_csv(files["flags"], sep="\t", dtype={"chrom": str}, low_memory=False)
    ann = pd.read_csv(files["annotation"], sep="\t", dtype={"chrom": str}, low_memory=False)
    assert panel["ModelID"].is_unique, "panel ModelIDs not unique"
    assert flags["gene_key"].is_unique and ann["depmap_column"].is_unique

    genes = flags.set_index("gene_key")
    a = ann.set_index("depmap_column")
    assert set(genes.index) == set(a.index), "frozen_flags and gene_annotation have different gene keys"
    a = a.loc[genes.index]
    assert (a["chrom"].fillna("NA").to_numpy() == genes["chrom"].fillna("NA").to_numpy()).all(), "chrom mismatch"
    genes["start"] = a["start"].to_numpy()
    genes["end"] = a["end"].to_numpy()
    genes["arm"] = genes["arm"].where(genes["arm"].isin(["p", "q"]))
    genes["arm_id"] = genes["chrom"].astype(str) + genes["arm"].astype(str)

    ids = panel["ModelID"].to_numpy()
    cn_df = pd.read_parquet(files["cn"]).set_index("ModelID")
    ex_df = pd.read_parquet(files["expr"]).set_index("ModelID")
    missing = set(ids) - set(cn_df.index) | set(ids) - set(ex_df.index)
    if missing:
        raise AssertionError(f"{len(missing)} panel models absent from matched matrices")

    # Context lines (G2/B4e) are taken from the full matched expression matrix, not the panel.
    context_model_check()
    context_expr = {}
    for name, (ach, _) in CONTEXT_LINES.items():
        if ach in ex_df.index:
            context_expr[name] = ex_df.loc[ach].astype(np.float64)

    cn_df = cn_df.loc[ids]
    ex_df = ex_df.loc[ids]
    # CN genes used by the engine: placed on 1-22/X, not PAR/Y/MT, finite in every panel line.
    in_scope = genes["chrom"].isin(CHROMS) & ~genes["is_par"] & ~genes["is_y"] & ~genes["is_mt"] & genes["arm"].notna()
    cn_cols = [c for c in cn_df.columns if c in genes.index and in_scope.get(c, False)]
    cn_block = cn_df[cn_cols]
    finite = np.isfinite(cn_block.to_numpy()).all(axis=0)
    cn_cols = list(np.asarray(cn_cols)[finite])
    cn_keys = gene_order(genes, cn_cols)
    relcn = cn_df[cn_keys].to_numpy(np.float64)

    expr_keys = np.asarray([c for c in ex_df.columns if c in genes.index])
    expr = ex_df[expr_keys].to_numpy(np.float64)
    for name in context_expr:
        context_expr[name] = context_expr[name].reindex(expr_keys)

    reg_mask = genes["variable_cn"] & in_scope
    tgt_mask = genes["measured"] & in_scope
    regulators = gene_order(genes, genes.index[reg_mask])
    targets = gene_order(genes, genes.index[tgt_mask])
    cn_set = set(cn_keys)
    if not set(regulators) <= cn_set:
        raise AssertionError("some variable-CN regulators have non-finite CN in the panel")
    if not set(targets) <= cn_set or not set(targets) <= set(expr_keys):
        raise AssertionError("some measured targets lack CN or expression")

    cks = {}
    if checksums:
        cks = {k: sha256(p) for k, p in files.items()}
    ds = Dataset(panel=panel.reset_index(drop=True), relcn=relcn, cn_keys=np.asarray(cn_keys), expr=expr,
                 expr_keys=expr_keys, genes=genes, regulators=np.asarray(regulators), targets=np.asarray(targets),
                 context_expr=context_expr, sample_set="S2", input_checksums=cks)
    return apply_sample_set(ds, sample_set)


CN_OUTLIER_FRAC = 0.01


def cn_outlier_score(ds: Dataset) -> np.ndarray:
    """Genome-wide max |log2 relCN| per line over the in-scope CN genes, with the design floor 2^-2 and no cap
    (DESIGN §13 CN-outlier sensitivity). The floor stops homozygous deletions (relCN ~ 1e-6 -> log2 ~ -20)
    from dominating; the score is then driven by high-level amplicons (e.g. COLO320 MYC)."""
    return np.abs(np.log2(np.maximum(ds.relcn, CN_FLOOR))).max(axis=1)


def apply_sample_set(ds: Dataset, sample_set: str) -> Dataset:
    """S2 (default), S1 (in_S1), S2_noHCT116, S2_noCNoutlier (top 1% of lines by cn_outlier_score removed)."""
    if sample_set in ("S2", "SIM") or sample_set == ds.sample_set:
        return ds
    if sample_set == "S1":
        return ds.subset_rows(ds.panel["in_S1"].astype(bool).to_numpy(), "S1")
    if sample_set == "S2_noHCT116":
        return ds.subset_rows(ds.panel["ModelID"].to_numpy() != CONTEXT_LINES["HCT116"][0], "S2_noHCT116")
    if sample_set == "S2_noCNoutlier":
        sc = cn_outlier_score(ds)
        k = int(np.ceil(CN_OUTLIER_FRAC * ds.n))
        drop = np.argsort(-sc, kind="stable")[:k]
        keep = np.ones(ds.n, bool)
        keep[drop] = False
        out = ds.subset_rows(keep, "S2_noCNoutlier")
        names = ds.panel["CellLineName"].to_numpy() if "CellLineName" in ds.panel else ds.model_ids
        out.sample_notes["cn_outlier_removed"] = [
            {"ModelID": str(ds.model_ids[i]), "name": str(names[i]), "max_abs_log2_relcn": float(sc[i])} for i in drop]
        out.sample_notes["cn_outlier_n_removed"] = int(k)
        out.sample_notes["cn_outlier_threshold"] = float(sc[drop].min())
        return out
    raise ValueError(f"unknown sample set {sample_set}")


# ----------------------------------------------------------------------------- invariants
def check_invariants(ds: Dataset, real: bool = True) -> dict:
    """DESIGN §8.6 invariants that concern the model inputs. Raises AssertionError on failure."""
    out = {}
    n = ds.n
    assert ds.relcn.shape == (n, len(ds.cn_keys)) and ds.expr.shape == (n, len(ds.expr_keys))
    assert len(ds.panel) == n and ds.panel["ModelID"].is_unique
    assert len(set(ds.cn_keys)) == len(ds.cn_keys) and len(set(ds.expr_keys)) == len(ds.expr_keys), "gene IDs not 1:1"
    assert np.isfinite(ds.relcn).all() and (ds.relcn > 0).all(), "CN has NaN/non-positive values"
    assert np.isfinite(ds.expr).all(), "expression has NaN"
    g = ds.genes
    assert not g.loc[ds.regulators, ["is_y", "is_par", "is_mt"]].any().any()
    assert not g.loc[ds.targets, ["is_y", "is_par", "is_mt"]].any().any()
    # median cis r positive (on the targets: x_B vs y_B)
    ci = pd.Index(ds.cn_keys).get_indexer(ds.targets)
    ei = pd.Index(ds.expr_keys).get_indexer(ds.targets)
    x = cn_transform(ds.relcn[:, ci])
    y = ds.expr[:, ei]
    x = x - x.mean(0)
    y = y - y.mean(0)
    r = (x * y).sum(0) / np.sqrt((x * x).sum(0) * (y * y).sum(0))
    out["median_cis_r"] = float(np.nanmedian(r))
    assert out["median_cis_r"] > 0, "median cis r is not positive"
    if real:
        sym = g.reset_index().set_index("symbol")
        for s, arm in {"ERBB2": "17q", "MYC": "8q", "CDKN2A": "9p"}.items():
            got = str(sym.loc[s, "chrom"]) + str(sym.loc[s, "arm"])
            assert got == arm, f"{s} on {got}, expected {arm}"
        for name, (ach, _) in CONTEXT_LINES.items():
            if name == "HCT116":
                assert name in ds.context_expr, "HCT116 context expression missing"
    out["n_samples"] = int(n)
    out["n_regulators"] = int(len(ds.regulators))
    out["n_targets"] = int(len(ds.targets))
    return out
