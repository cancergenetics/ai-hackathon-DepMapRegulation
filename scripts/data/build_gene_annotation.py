#!/usr/bin/env python
"""Gene annotation for DepMap gene columns and Perturb-seq (PIE X-Atlas/Orion) identifiers.

Inputs
  starting_data/gene_with_protein_product.txt      HGNC protein-coding set (shipped with project)
  data/external/annotation/hgnc_complete_set_<date>.tsv   HGNC complete set (all locus types)
  data/external/annotation/Homo_sapiens.GRCh38.<rel>.gtf.gz   Ensembl gene coordinates (GRCh38)
  data/external/annotation/cytoBand.hg38.txt.gz    UCSC cytobands (arm / band by gene midpoint)
  data/interim/depmap_{cn_wgs,expr_tpm_logp1}.parquet   (column names only)
  data/processed/pie/{vocab_targets,vocab_features,targets_*,features_*}.tsv

Outputs
  data/processed/gene_annotation.tsv   one row per DepMap gene column (union of CN and expression)
  data/processed/pie_gene_map.tsv      one row per PIE identifier (union of targets and features)
  data/processed/gene_annotation_summary.json   counts used in docs/DATA_PROFILE.md

Mapping rules (explicit; see docs/DATA_PROFILE.md):
  DepMap "SYMBOL (ENTREZ)": Entrez -> HGNC (complete set). If the Entrez id is not in HGNC, fall back to
    the symbol (current > previous > alias; only unique hits accepted).
  DepMap bare "ENSG..." columns (CN only): Ensembl id -> HGNC; "_PAR_Y" suffix = Y-chromosome copy of a
    pseudoautosomal gene (coordinates translated from X using GRCh38 PAR offsets).
  PIE identifiers: "ENSG..." -> Ensembl id; "SYMBOL_ENSG..." (10x duplicate-name convention) -> embedded
    Ensembl id; symbols -> HGNC current symbol > Ensembl gene_name (unique) >
    HGNC previous symbol (unique) > HGNC alias (unique); otherwise unmapped / ambiguous.
  Coordinates/chromosome come from the Ensembl GTF (via the HGNC Ensembl id, else a unique GTF gene_name
    match). Arm (p/q) and band come from the UCSC band containing the gene midpoint. Genes without
    a primary-assembly location fall back to the HGNC 'location' field for chromosome/arm.
  chrom_class: autosome / X / Y / MT / unplaced (scaffold or no reference-assembly location).
  PAR genes: is_par (HGNC location lists X and Y, or coordinates in a GRCh38 PAR), par_region
    (PAR1 / PAR2 / hgnc_location_only), par_y_copy (DepMap "_PAR_Y" CN column). In pie_gene_map.tsv,
    depmap_cn_column is always the X/autosomal copy; the Y copy is in depmap_cn_column_par_y.
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
ANN = ROOT / "data" / "external" / "annotation"
PROC = ROOT / "data" / "processed"
PIE = PROC / "pie"
HGNC_FILE = ANN / "hgnc_complete_set_2026-10-06.tsv"
GTF_FILE = ANN / "Homo_sapiens.GRCh38.116.gtf.gz"
CYTO_FILE = ANN / "cytoBand.hg38.txt.gz"
META = ["SequencingID", "ModelConditionID", "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel"]
PRIMARY = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
# GRCh38 pseudoautosomal regions (X coords -> Y coords). PAR1 has identical coordinates on X and Y.
PAR = [(10_001, 2_781_479, 10_001), (155_701_383, 156_030_895, 56_887_903)]
SYM_ENSG = re.compile(r"^.+_(ENSG\d{11})$")
DEPMAP_COL = re.compile(r"^(?P<sym>.+) \((?P<entrez>\d+)\)$")
ENSG = re.compile(r"^(ENSG\d{11})(?:\.\d+)?(?P<par>_PAR_Y)?$")


# ----------------------------------------------------------------------------------------------
def load_hgnc() -> pd.DataFrame:
    h = pd.read_csv(HGNC_FILE, sep="\t", dtype=str, low_memory=False)
    h = h[h.status == "Approved"].copy()
    pc = pd.read_csv(ROOT / "starting_data" / "gene_with_protein_product.txt", sep="\t", dtype=str)
    missing = set(pc.hgnc_id) - set(h.hgnc_id)
    if missing:
        print(f"NOTE: {len(missing)} protein-coding HGNC ids from starting_data not in complete set "
              f"(newer file); appending them", file=sys.stderr)
        h = pd.concat([h, pc[pc.hgnc_id.isin(missing)]], ignore_index=True)
    return h


def split_multi(s: str | float) -> list[str]:
    if not isinstance(s, str) or not s:
        return []
    return [x.strip() for x in s.split("|") if x.strip()]


def load_gtf() -> pd.DataFrame:
    rows = []
    attr = re.compile(r'(\w+) "([^"]*)"')
    with gzip.open(GTF_FILE, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if parts[2] != "gene":
                continue
            a = dict(attr.findall(parts[8]))
            rows.append((a["gene_id"], a.get("gene_name", ""), a.get("gene_biotype", ""), parts[0],
                         int(parts[3]), int(parts[4]), parts[6]))
    g = pd.DataFrame(rows, columns=["ensembl_gene_id", "gtf_gene_name", "gtf_biotype", "gtf_seqname",
                                    "start", "end", "strand"])
    return g


def load_cyto() -> pd.DataFrame:
    c = pd.read_csv(CYTO_FILE, sep="\t", header=None, names=["chrom", "s", "e", "band", "stain"])
    c["chrom"] = c.chrom.str.replace("^chr", "", regex=True)
    c = c[c.chrom.isin(PRIMARY)]
    return c


def band_at(cyto: pd.DataFrame, chrom: str, pos: int) -> str | None:
    sub = cyto[cyto.chrom == chrom]
    hit = sub[(sub.s < pos) & (sub.e >= pos)]
    return None if hit.empty else hit.band.iloc[0]


def hgnc_loc_parse(loc: str | float) -> tuple[str | None, str | None]:
    """Chromosome and arm from an HGNC location string such as '17q12', 'Xp22.33 and Yp11.2'."""
    if not isinstance(loc, str):
        return None, None
    if loc.startswith("mitochondria"):
        return "MT", None
    m = re.match(r"^(\d{1,2}|X|Y)(?:\s*(?:([pq])|cen))?", loc)
    if not m:
        return None, None
    return m.group(1), m.group(2)


# ----------------------------------------------------------------------------------------------
class Mapper:
    def __init__(self, hgnc: pd.DataFrame, gtf: pd.DataFrame):
        self.h = hgnc.set_index("hgnc_id")
        self.by_entrez = {e: i for i, e in zip(hgnc.hgnc_id, hgnc.entrez_id) if isinstance(e, str)}
        self.by_ens = {e: i for i, e in zip(hgnc.hgnc_id, hgnc.ensembl_gene_id) if isinstance(e, str)}
        self.by_sym = {s: i for i, s in zip(hgnc.hgnc_id, hgnc.symbol)}
        self.by_prev = defaultdict(set)
        self.by_alias = defaultdict(set)
        for i, p, a in zip(hgnc.hgnc_id, hgnc.prev_symbol, hgnc.alias_symbol):
            for s in split_multi(p):
                self.by_prev[s].add(i)
            for s in split_multi(a):
                self.by_alias[s].add(i)
        self.gtf = gtf.set_index("ensembl_gene_id")
        prim = gtf[gtf.gtf_seqname.isin(PRIMARY)]
        cnt = prim.gtf_gene_name.value_counts()
        uniq = prim[prim.gtf_gene_name.isin(cnt[cnt == 1].index) & (prim.gtf_gene_name != "")]
        self.gtf_by_name = dict(zip(uniq.gtf_gene_name, uniq.ensembl_gene_id))

    def by_symbol(self, sym: str) -> tuple[str | None, str, str]:
        """Return (hgnc_id, method, note)."""
        if sym in self.by_sym:
            note = ""
            if sym in self.by_prev or sym in self.by_alias:
                note = "also prev/alias of other gene(s)"
            return self.by_sym[sym], "hgnc_symbol", note
        if sym in self.gtf_by_name:
            ens = self.gtf_by_name[sym]
            if ens in self.by_ens:
                return self.by_ens[ens], "ensembl_gene_name", ""
            return None, "ensembl_gene_name_no_hgnc", ens
        prev = self.by_prev.get(sym, set())
        if len(prev) == 1:
            return next(iter(prev)), "hgnc_prev_symbol", ""
        if len(prev) > 1:
            return None, "ambiguous_prev_symbol", ";".join(sorted(prev))
        alias = self.by_alias.get(sym, set())
        if len(alias) == 1:
            return next(iter(alias)), "hgnc_alias_symbol", ""
        if len(alias) > 1:
            return None, "ambiguous_alias_symbol", ";".join(sorted(alias))
        return None, "unmapped", ""


def locate(m: Mapper, cyto: pd.DataFrame, hgnc_id: str | None, ens: str | None, par_y: bool = False) -> dict:
    out = dict(hgnc_id=hgnc_id, symbol=None, entrez_id=None, ensembl_gene_id=ens, locus_group=None,
               locus_type=None, hgnc_location=None, gtf_gene_name=None, gtf_biotype=None, chrom=None,
               start=np.nan, end=np.nan, strand=None, arm=None, cytoband=None, chrom_class=None,
               coord_source=None, is_par=False)
    if hgnc_id is not None:
        r = m.h.loc[hgnc_id]
        out.update(symbol=r.symbol, entrez_id=r.entrez_id, locus_group=r.locus_group, locus_type=r.locus_type,
                   hgnc_location=r.location)
        if ens is None and isinstance(r.ensembl_gene_id, str):
            ens = r.ensembl_gene_id
            out["ensembl_gene_id"] = ens
        if isinstance(r.location, str) and "X" in r.location and "Y" in r.location and "and" in r.location:
            out["is_par"] = True
    if ens is None and out["symbol"] is not None and out["symbol"] in m.gtf_by_name:
        ens = m.gtf_by_name[out["symbol"]]
        out["ensembl_gene_id"] = ens
    if ens is not None and ens in m.gtf.index:
        g = m.gtf.loc[ens]
        out.update(gtf_gene_name=g.gtf_gene_name, gtf_biotype=g.gtf_biotype)
        if g.gtf_seqname in PRIMARY:
            chrom, s, e = g.gtf_seqname, int(g.start), int(g.end)
            if par_y:
                for xs, xe, ys in PAR:
                    if xs <= s and e <= xe:
                        chrom, s, e = "Y", s - xs + ys, e - xs + ys
                        break
                else:
                    chrom = None
            if chrom is not None:
                out.update(chrom=chrom, start=s, end=e, strand=g.strand, coord_source=f"ensembl_{GTF_FILE.name.split('.')[2]}")
                band = band_at(cyto, chrom, (s + e) // 2) if chrom != "MT" else None
                out.update(cytoband=(chrom + band) if band else None, arm=band[0] if band else None)
                if any(xs <= s and e <= xe for xs, xe, _ in PAR) and chrom == "X":
                    out["is_par"] = True
                if par_y:
                    out["is_par"] = True
        else:
            out["coord_source"] = f"scaffold:{g.gtf_seqname}"
    if out["chrom"] is None:
        c, a = hgnc_loc_parse(out["hgnc_location"])
        if par_y and c == "X":
            c = "Y"
        if c is not None and out["coord_source"] is None:
            out.update(chrom=c, arm=a, coord_source="hgnc_location")
    c = out["chrom"]
    if c is None or (isinstance(out["coord_source"], str) and out["coord_source"].startswith("scaffold")):
        out["chrom_class"] = "unplaced"
        out["chrom"] = None if out["coord_source"] is None or out["coord_source"].startswith("scaffold") else c
    elif c in ("X", "Y", "MT"):
        out["chrom_class"] = c
    else:
        out["chrom_class"] = "autosome"
    return out


def par_region(chrom, start, end) -> str:
    """PAR1/PAR2 from GRCh38 coordinates; 'hgnc_location_only' when HGNC lists X and Y but the gene
    lies outside the PARs (X-transposed / X-degenerate homolog pairs)."""
    if chrom in ("X", "Y") and pd.notna(start):
        for (xs, xe, ys), name in zip(PAR, ("PAR1", "PAR2")):
            lo, hi = (xs, xe) if chrom == "X" else (ys, ys + (xe - xs))
            if lo <= start and end <= hi:
                return name
    return "hgnc_location_only"


# ----------------------------------------------------------------------------------------------
def depmap_annotation(m: Mapper, cyto: pd.DataFrame) -> pd.DataFrame:
    cn = [c for c in pl.read_parquet_schema(ROOT / "data/interim/depmap_cn_wgs.parquet") if c not in META]
    ex = [c for c in pl.read_parquet_schema(ROOT / "data/interim/depmap_expr_tpm_logp1.parquet") if c not in META]
    cols = sorted(set(cn) | set(ex))
    cn_s, ex_s = set(cn), set(ex)
    rows = []
    for col in cols:
        rec = dict(depmap_column=col, in_cn=col in cn_s, in_expr=col in ex_s, depmap_symbol=None,
                   depmap_entrez=None, map_method=None, map_note="")
        hid, ens, par_y = None, None, False
        mm = DEPMAP_COL.match(col)
        me = ENSG.match(col)
        if mm:
            sym, ent = mm.group("sym"), mm.group("entrez")
            rec.update(depmap_symbol=sym, depmap_entrez=ent)
            if ent in m.by_entrez:
                hid = m.by_entrez[ent]
                rec["map_method"] = "entrez"
                cur = m.h.loc[hid, "symbol"]
                if cur != sym:
                    rec["map_note"] = f"symbol updated {sym}->{cur}"
            else:
                hid, meth, note = m.by_symbol(sym)
                rec.update(map_method=f"entrez_not_in_hgnc;{meth}", map_note=note)
        elif me:
            ens, par_y = me.group(1), bool(me.group("par"))
            hid = m.by_ens.get(ens)
            rec["map_method"] = "ensembl" + ("_PAR_Y" if par_y else "")
            if hid is None:
                rec["map_note"] = "ensembl id not in HGNC"
        else:
            rec["map_method"] = "unparsed"
        rec.update(locate(m, cyto, hid, ens, par_y))
        rows.append(rec)
    df = pd.DataFrame(rows)
    df["par_y_copy"] = df.depmap_column.str.endswith("_PAR_Y")
    df["par_region"] = [par_region(r.chrom, r.start, r.end) if r.is_par else None for r in df.itertuples()]
    dup = df.hgnc_id.notna() & df.duplicated(["hgnc_id", "chrom"], keep=False)
    df["hgnc_dup_in_depmap"] = dup
    return df


def pie_map(m: Mapper, cyto: pd.DataFrame, dep: pd.DataFrame) -> pd.DataFrame:
    tv = pd.read_csv(PIE / "vocab_targets.tsv", sep="\t", dtype={"target": str}, keep_default_na=False)
    fv = pd.read_csv(PIE / "vocab_features.tsv", sep="\t", dtype={"feature": str}, keep_default_na=False)
    ids = sorted(set(tv.target) | set(fv.feature))
    roles = {}
    for cell in ("hct116", "hek293t"):
        t = pd.read_csv(PIE / f"targets_{cell}.tsv", sep="\t", dtype={"target": str}, keep_default_na=False)
        f = pd.read_csv(PIE / f"features_{cell}.tsv", sep="\t", dtype={"feature": str}, keep_default_na=False)
        roles[f"target_{cell}"] = set(t.target)
        roles[f"feature_{cell}"] = set(f.feature)
    # DepMap lookup by hgnc id (primary key) -> columns
    # One DepMap column per HGNC id. The Y copies of PAR genes ("ENSG..._PAR_Y", CN only) are linked in a
    # separate column so that depmap_cn_column always refers to the X/autosomal copy (data review E4).
    keep = dep.hgnc_id.notna() & ~dep.par_y_copy
    dep_cn = dep[dep.in_cn & keep].groupby("hgnc_id").depmap_column.agg(";".join)
    dep_ex = dep[dep.in_expr & keep].groupby("hgnc_id").depmap_column.agg(";".join)
    dep_cn_y = dep[dep.in_cn & dep.hgnc_id.notna() & dep.par_y_copy].groupby("hgnc_id").depmap_column.agg(";".join)
    for name, ser in [("cn", dep_cn), ("expr", dep_ex), ("cn_par_y", dep_cn_y)]:
        if ser.str.contains(";").any():
            raise SystemExit(f"multiple DepMap {name} columns for one HGNC id: {ser[ser.str.contains(';')].to_dict()}")
    rows = []
    for ident in ids:
        rec = dict(pie_id=ident)
        for k, s in roles.items():
            rec[k] = ident in s
        me = ENSG.match(ident)
        msx = SYM_ENSG.match(ident)
        if msx and msx.group(1) in m.gtf.index:
            # 10x-style "SYMBOL_ENSG..." name used when one symbol has several gene ids
            ens = msx.group(1)
            hid = m.by_ens.get(ens)
            rec.update(id_type="symbol_ensembl", map_method="embedded_ensembl", map_note="" if hid else "no HGNC")
            loc = locate(m, cyto, hid, ens)
        elif me:
            ens = me.group(1)
            hid = m.by_ens.get(ens)
            rec.update(id_type="ensembl", map_method="ensembl" if ens in m.gtf.index else "ensembl_not_in_gtf",
                       map_note="" if hid else "no HGNC")
            loc = locate(m, cyto, hid, ens)
        else:
            hid, meth, note = m.by_symbol(ident)
            rec.update(id_type="symbol", map_method=meth, map_note=note)
            ens = note if meth == "ensembl_gene_name_no_hgnc" else None
            loc = locate(m, cyto, hid, ens)
        rec.update(loc)
        rec["depmap_cn_column"] = dep_cn.get(hid) if hid else None
        rec["depmap_expr_column"] = dep_ex.get(hid) if hid else None
        rec["depmap_cn_column_par_y"] = dep_cn_y.get(hid) if hid else None
        rows.append(rec)
    df = pd.DataFrame(rows)
    # one-to-one: if several PIE ids map to the same HGNC id (typically HGNC locus merges, where a
    # previous symbol is still a separate feature in the PIE reference), keep the highest-priority id
    # and mark the others as collisions (location kept, DepMap link removed).
    prio = {"hgnc_symbol": 0, "ensembl": 1, "embedded_ensembl": 1, "ensembl_gene_name": 2, "hgnc_prev_symbol": 3, "hgnc_alias_symbol": 4}
    df["collision_with"] = None
    df["_p"] = df.map_method.map(prio).fillna(9)
    for hid, grp in df[df.hgnc_id.notna()].groupby("hgnc_id"):
        if len(grp) < 2:
            continue
        grp = grp.sort_values(["_p", "pie_id"])
        keep = grp.index[0]
        for i in grp.index[1:]:
            df.loc[i, "collision_with"] = df.loc[keep, "pie_id"]
            df.loc[i, "map_note"] = (str(df.loc[i, "map_note"]) + f" collision:{df.loc[i, 'map_method']}").strip()
            df.loc[i, ["depmap_cn_column", "depmap_expr_column", "depmap_cn_column_par_y"]] = None
            df.loc[i, "map_method"] = "collision"
    return df.drop(columns="_p")


def summarise(dep: pd.DataFrame, pie: pd.DataFrame) -> dict:
    s: dict = {}
    s["depmap"] = {
        "n_columns": int(len(dep)), "n_cn": int(dep.in_cn.sum()), "n_expr": int(dep.in_expr.sum()),
        "n_both": int((dep.in_cn & dep.in_expr).sum()),
        "map_method": dep.map_method.value_counts().to_dict(),
        "n_symbol_updated": int(dep.map_note.str.startswith("symbol updated").sum()),
        "n_no_hgnc": int(dep.hgnc_id.isna().sum()),
        "chrom_class": dep.chrom_class.value_counts().to_dict(),
        "chrom_class_cn": dep[dep.in_cn].chrom_class.value_counts().to_dict(),
        "chrom_class_expr": dep[dep.in_expr].chrom_class.value_counts().to_dict(),
        "coord_source": dep.coord_source.fillna("none").value_counts().to_dict(),
        "n_with_coords": int(dep.start.notna().sum()),
        "n_arm_missing": int(dep.arm.isna().sum()),
        "n_par": int(dep.is_par.sum()),
        "par_region": dep.par_region.dropna().value_counts().to_dict(),
        "par_genes": sorted(dep[dep.is_par].depmap_column.tolist()),
        "n_hgnc_dup": int(dep.hgnc_dup_in_depmap.sum()),
        "locus_group": dep.locus_group.fillna("none").value_counts().to_dict(),
    }
    # chromosome agreement HGNC location vs Ensembl
    hl = dep.hgnc_location.map(lambda x: hgnc_loc_parse(x)[0])
    both = dep.chrom.notna() & hl.notna() & dep.coord_source.astype(str).str.startswith("ensembl")
    s["depmap"]["hgnc_vs_ensembl_chrom_disagree"] = int((both & (hl != dep.chrom) & ~dep.is_par).sum())
    ha = dep.hgnc_location.map(lambda x: hgnc_loc_parse(x)[1])
    both_a = both & ha.notna() & dep.arm.notna() & (hl == dep.chrom)
    s["depmap"]["hgnc_vs_ucsc_arm_disagree"] = int((both_a & (ha != dep.arm)).sum())
    out = {}
    for role in ["target_hct116", "target_hek293t", "feature_hct116", "feature_hek293t"]:
        sub = pie[pie[role]]
        out[role] = {
            "n": int(len(sub)),
            "id_type": sub.id_type.value_counts().to_dict(),
            "map_method": sub.map_method.value_counts().to_dict(),
            "n_mapped_hgnc": int(sub.hgnc_id.notna().sum()),
            "n_with_chrom": int(sub.chrom.notna().sum()),
            "chrom_class": sub.chrom_class.value_counts().to_dict(),
            "n_in_depmap_cn": int(sub.depmap_cn_column.notna().sum()),
            "n_in_depmap_expr": int(sub.depmap_expr_column.notna().sum()),
            "n_in_depmap_both": int((sub.depmap_cn_column.notna() & sub.depmap_expr_column.notna()).sum()),
            "locus_group": sub.locus_group.fillna("none").value_counts().to_dict(),
        }
    s["pie"] = out
    dup = pie[pie.hgnc_id.notna()].groupby("hgnc_id").pie_id.agg(list)
    s["pie_n_collisions"] = int((pie.map_method == "collision").sum())
    s["pie_rescues"] = {m: sorted(pie[pie.map_method == m].pie_id.tolist())
                        for m in ["hgnc_prev_symbol", "hgnc_alias_symbol", "ensembl_gene_name", "ensembl_gene_name_no_hgnc",
                                  "ambiguous_prev_symbol", "ambiguous_alias_symbol", "unmapped", "collision",
                                  "embedded_ensembl", "ensembl_not_in_gtf"]}
    s["pie_many_to_one_hgnc"] = {k: v for k, v in dup.items() if len(v) > 1}
    return s


def main() -> int:
    hgnc = load_hgnc()
    gtf = load_gtf()
    cyto = load_cyto()
    m = Mapper(hgnc, gtf)
    dep = depmap_annotation(m, cyto)
    dep.to_csv(PROC / "gene_annotation.tsv", sep="\t", index=False)
    pie = pie_map(m, cyto, dep)
    pie.to_csv(PROC / "pie_gene_map.tsv", sep="\t", index=False)
    s = summarise(dep, pie)
    s["sources"] = {"hgnc": HGNC_FILE.name, "gtf": GTF_FILE.name, "cytoband": CYTO_FILE.name}
    with open(PROC / "gene_annotation_summary.json", "w") as fh:
        json.dump(s, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in s.items() if k != "pie_many_to_one_hgnc"}, indent=1, default=str)[:6000])
    print("many-to-one PIE ids:", len(s["pie_many_to_one_hgnc"]), list(s["pie_many_to_one_hgnc"].items())[:10])
    return 0


if __name__ == "__main__":
    sys.exit(main())
