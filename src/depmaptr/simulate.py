"""SIM (DESIGN §7): toy data with known truth, in the same `Dataset` structure as the real data.

Genome: autosomes 1..n_auto (arms p, q of `genes_per_arm` genes in segments of 3-8 genes) and X (two arms).
CN (log2 relative):  arm event + focal segment event + gene noise (sd 0.05); X: male -1.
Each confounder below has a test in tests/test_sim.py that passes with the adjustment and FAILS when the
adjustment is removed (ablation tests):
- **cross-chromosome arm co-selection** (needs the LOCO arm block): `coselect` arms (1q, 3p) move together,
  up or down with equal magnitude (so |arm|-based burden cannot absorb it), and 3p dosage drives
  `driver_targets` (an arm-level driver gene on 3p);
- **lineage confounding** (needs lineage in Z): lineage L1 has frequent `lineage_arm` (2q) gains AND an
  expression programme on `lineage_programme` genes that is not caused by CN;
- **sex confounding** (needs sex in Z): raw X CN is -1 in males, and an autosomal `sex_programme` is up in
  males; one planted regulator is on X;
- **homozygous deletions** (needs the deletion indicator in W_B): segments in `homdel_arm` are deleted
  (relCN = 1e-6) in ~4% of lines, with collapsed expression (non-linear in x); a partner focal segment on
  `homdel_partner_chrom`q is amplified in most of those lines;
- **co-segment neighbours** (A3): planted regulators sit in focal segments whose other genes share their CN
  but have weak cis expression; half of the planted segments contain a CN-only neighbour;
- **negative-cis regulator** (needs the s_A rule): cis slope -1, trans effects through its mRNA;
- **CN-only genes**: `cn_only_frac` of genes have no expression column (not targets; t_cis = 0).

Truth: `truth["planted"]` DataFrame (regulator, target, gamma, kd_sign), where kd_sign is the sign of the
change in y_B when A's mRNA is lowered (CRISPRi).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .io import Dataset, gene_order

LINEAGES = ["L1", "L2", "L3", "L4", "L5", "L6"]


def simulate(seed: int = 1, n: int = 600, n_auto: int = 6, genes_per_arm: int = 60, x_genes_per_arm: int = 25,
             null: bool = False, n_planted: int = 12, targets_per_reg: int = 25, gamma: float = 0.5,
             coselect=("1q", "3p"), lineage_arm: str = "2q", driver_arm: str = "3p",
             homdel_arm: str = "4p", homdel_partner_chrom: str = "5", noise: float = 0.5,
             n_lineage_programme: int = 40, n_sex_programme: int = 30, sex_programme_chrom: str = "6",
             driver_effect: float = 1.0, cn_only_frac: float = 0.06) -> tuple[Dataset, dict]:
    rng = np.random.default_rng(seed)
    chroms = [str(i) for i in range(1, n_auto + 1)] + ["X"]
    # ------------------------------------------------------------------ genome
    rows = []
    seg_id = 0
    for c in chroms:
        gpa = x_genes_per_arm if c == "X" else genes_per_arm
        pos = 0
        for a in ("p", "q"):
            k = 0
            while k < gpa:
                L = int(min(rng.integers(3, 9), gpa - k))
                for j in range(L):
                    rows.append((f"G{c}{a}{k + j:03d} ({len(rows) + 1})", c, a, pos, seg_id))
                    pos += 1000
                k += L
                seg_id += 1
            pos += 10_000_000
    genes = pd.DataFrame(rows, columns=["gene_key", "chrom", "arm", "start", "segment"]).set_index("gene_key")
    genes["arm_id"] = genes["chrom"] + genes["arm"]
    G = len(genes)
    n_seg = seg_id
    keys = genes.index.to_numpy()
    arm_ids = sorted(genes["arm_id"].unique())
    seg_of = genes["segment"].to_numpy()
    seg_arm = genes.groupby("segment")["arm_id"].first()
    seg_size = genes.groupby("segment").size()

    # ------------------------------------------------------------------ samples
    lineage = rng.choice(LINEAGES, size=n, p=[0.2, 0.2, 0.15, 0.15, 0.15, 0.15])
    sex = rng.choice(["Male", "Female"], size=n)
    male = sex == "Male"

    # arm events. Co-selected arms move TOGETHER in both directions with equal magnitude (+-0.58), so the
    # genome-burden covariates (|arm|-based) cannot absorb the co-selection; only the LOCO arm block can.
    arm_lv = pd.DataFrame(0.0, index=range(n), columns=arm_ids)
    latent = rng.choice([1.0, -1.0, 0.0], size=n, p=[0.25, 0.25, 0.5])
    for aid in arm_ids:
        u = rng.random(n)
        if aid in coselect:
            arm_lv[aid] = np.where((latent != 0) & (u < 0.85), 0.58 * latent, 0.0)
            continue
        pg = np.full(n, 0.10)
        pl = np.full(n, 0.10)
        if aid == lineage_arm:
            pg = np.where(lineage == "L1", 0.7, 0.08)
        arm_lv[aid] = np.where(u < pg, 0.58, np.where(u < pg + pl, -1.0, 0.0))
    # focal segment events
    seg_lv = np.zeros((n, n_seg))
    u = rng.random((n, n_seg))
    seg_lv[u < 0.04] = rng.uniform(1.0, 2.0, size=(u < 0.04).sum())
    # focal loss -0.9 (not -1): arm loss + focal loss then stays above the 0.25 homdel threshold, so the
    # deletion indicator marks only true homozygous deletions in SIM
    seg_lv[(u >= 0.04) & (u < 0.08)] = -0.9
    # homozygous deletions: 2 segments of homdel_arm deleted in ~4% of lines; a partner focal segment on
    # another chromosome is amplified in most of those lines (cross-chromosome co-selection with the deletion).
    hd_segs = genes.loc[genes["arm_id"] == homdel_arm, "segment"].unique()[1:3]
    hd_lines = rng.random(n) < 0.04
    partner_seg = [s for s in seg_size.index if seg_arm[s] == homdel_partner_chrom + "q" and seg_size[s] >= 4][0]
    seg_lv[hd_lines & (rng.random(n) < 0.8), partner_seg] = 1.5
    log2cn = arm_lv[genes["arm_id"].to_numpy()].to_numpy() + seg_lv[:, seg_of]
    log2cn += rng.normal(0, 0.05, size=(n, G))
    isx = (genes["chrom"] == "X").to_numpy()
    log2cn[np.ix_(male, isx)] -= 1.0
    relcn = 2.0 ** log2cn
    hd_mask = np.zeros((n, G), bool)
    hd_mask[np.ix_(hd_lines, np.isin(seg_of, hd_segs))] = True
    relcn[hd_mask] = 1e-6
    x = np.log2(np.clip(relcn, 0.25, 8.0))

    # ------------------------------------------------------------------ expression
    mu = rng.uniform(2, 6, G)
    b_cis = rng.uniform(0.0, 0.4, G)                     # weak default cis dosage
    lin_eff = np.zeros((n, G))
    lp_genes = rng.choice(np.flatnonzero(genes["chrom"].isin(["4", "5"]).to_numpy() & ~np.isin(seg_of, hd_segs)
                                         & (seg_of != partner_seg)), size=n_lineage_programme, replace=False)
    lin_eff[np.ix_(lineage == "L1", lp_genes)] += 1.5
    for L in LINEAGES[1:]:                              # mild lineage shifts on all genes
        lin_eff[lineage == L] += rng.normal(0, 0.2, G)
    # sex programme: autosomal genes up in males (a sex confounder for X regulators, whose raw CN tracks sex)
    sp_genes = rng.choice(np.flatnonzero((genes["chrom"] == sex_programme_chrom).to_numpy()), size=n_sex_programme,
                          replace=False)
    lin_eff[np.ix_(male, sp_genes)] += 1.0
    xc = x - np.where(isx[None, :] & male[:, None], -1.0, 0.0)   # sex-centred dosage drives expression
    e = mu + lin_eff + rng.normal(0, noise, size=(n, G))

    # planted regulators: one gene in each of several multi-gene segments (one on X), strong cis.
    banned_arms = set(coselect) | {lineage_arm, homdel_arm}
    eligible = [s for s in seg_size.index if seg_size[s] >= 4 and seg_arm[s] not in banned_arms
                and not seg_arm[s].startswith("X") and s != partner_seg and seg_arm[s][:-1] != sex_programme_chrom]
    x_eligible = [s for s in seg_size.index if seg_size[s] >= 4 and seg_arm[s].startswith("X")]
    planted_segs = list(rng.choice(eligible, size=n_planted - 1, replace=False)) + [rng.choice(x_eligible)]
    planted_regs = np.array([rng.choice(np.flatnonzero(seg_of == s)) for s in planted_segs])
    neg_cis_reg = planted_regs[0]
    x_planted = planted_regs[-1]
    b_cis[planted_regs] = 1.0
    b_cis[neg_cis_reg] = -1.0
    for s in planted_segs:                              # weak-cis co-segment neighbours
        for g in np.flatnonzero(seg_of == s):
            if g not in planted_regs:
                b_cis[g] = 0.1
    driver = np.flatnonzero((genes["arm_id"] == driver_arm).to_numpy())[5]
    y = e + b_cis * xc
    y = np.where(hd_mask, mu * 0.1, y)                  # homozygous deletion: expression collapses

    truth_rows = []
    excluded_targets = (set(lp_genes) | set(sp_genes) | set(planted_regs) | {driver}
                        | set(np.flatnonzero(np.isin(seg_of, hd_segs)))
                        | set(np.flatnonzero(np.isin(seg_of, planted_segs))))
    if not null:
        used = set()
        for A in planted_regs:
            cand = np.flatnonzero((genes["chrom"] != genes["chrom"].iloc[A]).to_numpy())
            cand = np.array([b for b in cand if b not in excluded_targets and b not in used and not isx[b]])
            tg = rng.choice(cand, size=targets_per_reg, replace=False)
            used |= set(tg)
            for b in tg:
                gm = gamma * rng.choice([-1.0, 1.0]) * rng.uniform(0.6, 1.4)
                y[:, b] += gm * (y[:, A] - y[:, A].mean())
                truth_rows.append((keys[A], keys[b], gm, -np.sign(gm)))
        # arm-level driver effects (3p dosage -> targets elsewhere); not in `planted` (arm-level truth)
        cand = np.flatnonzero(~genes["chrom"].isin([driver_arm[:-1], "X"]).to_numpy())
        cand = np.array([b for b in cand if b not in excluded_targets and b not in used])
        driver_targets = rng.choice(cand, size=20, replace=False)
        for b in driver_targets:
            y[:, b] += driver_effect * xc[:, driver]
        used |= set(driver_targets)
    else:
        driver_targets = np.array([], int)
        used = set()
    expr = y

    # CN-only genes (no expression column): random genes, plus one weak neighbour in half of the planted
    # segments; never planted regulators, their targets, or confounder genes.
    protected = set(planted_regs) | used | excluded_targets
    pool = np.array([g for g in range(G) if g not in protected])
    cn_only = set(rng.choice(pool, size=int(cn_only_frac * G), replace=False))
    for s in planted_segs[: len(planted_segs) // 2]:
        nb = [g for g in np.flatnonzero(seg_of == s) if g not in planted_regs and g not in used]
        if nb:
            cn_only.add(nb[0])
    has_expr = np.ones(G, bool)
    has_expr[list(cn_only)] = False

    panel = pd.DataFrame({
        "ModelID": [f"SIM-{i:04d}" for i in range(n)], "lineage": lineage, "primary_disease": lineage,
        "sex": sex, "ploidy": rng.normal(3, 0.7, n), "ploidy_missing": False, "wgd": (rng.random(n) < 0.6) * 1.0,
        "wgd_missing": False, "msi_score": rng.gamma(1.0, 3.0, n), "msi_missing": False,
        "tp53_altered": rng.random(n) < 0.6, "tp53_missing": False, "media_base": rng.choice(["RPMI", "DMEM"], n),
        "serum_free": False, "growth_pattern": rng.choice(["Adherent", "Suspension"], n, p=[0.8, 0.2]),
        "stranded": rng.random(n) < 0.5, "in_S1": True})
    miss = rng.random(n) < 0.2                          # missingness to exercise indicators
    panel.loc[miss, ["ploidy", "wgd"]] = np.nan
    panel.loc[miss, ["ploidy_missing", "wgd_missing"]] = True

    # flags: variable CN (>= 20 lines |log2| >= 1; X sex-centred), measured (CN and expression, SD > 0)
    lc = np.log2(relcn)
    lc_c = lc.copy()
    for sx in ("Male", "Female"):
        m = sex == sx
        lc_c[np.ix_(m, isx)] -= np.median(lc[np.ix_(m, isx)], axis=0)
    nvar = (np.abs(np.where(isx[None, :], lc_c, lc)) >= 1).sum(0)
    genes["variable_cn"] = nvar >= 20
    genes["has_cn"] = True
    genes["has_expr"] = has_expr
    genes["measured"] = has_expr & (expr.std(0) > 0)
    for f in ("is_par", "is_y", "is_mt"):
        genes[f] = False
    genes["symbol"] = [k.split(" ")[0] for k in keys]
    order = gene_order(genes, keys)
    pos = {k: i for i, k in enumerate(keys)}
    oi = np.array([pos[k] for k in order])
    eorder = np.array([k for k in order if has_expr[pos[k]]])
    ei = np.array([pos[k] for k in eorder])
    ds = Dataset(panel=panel, relcn=relcn[:, oi], cn_keys=order, expr=expr[:, ei], expr_keys=eorder, genes=genes,
                 regulators=gene_order(genes, genes.index[genes["variable_cn"]]),
                 targets=gene_order(genes, genes.index[genes["measured"]]),
                 context_expr={"HCT116": pd.Series(expr[0, ei], index=eorder)}, sample_set="SIM")
    truth = {
        "planted": pd.DataFrame(truth_rows, columns=["regulator", "target", "gamma", "kd_sign"]),
        "planted_regulators": keys[planted_regs], "neg_cis_regulator": keys[neg_cis_reg],
        "x_planted_regulator": keys[x_planted],
        "planted_segments": {keys[A]: list(keys[seg_of == s]) for A, s in zip(planted_regs, planted_segs)},
        "driver": keys[driver], "driver_targets": keys[driver_targets], "driver_arm": driver_arm,
        "coselect": coselect, "lineage_arm": lineage_arm, "lineage_programme": keys[lp_genes],
        "sex_programme": keys[sp_genes], "homdel_genes": keys[hd_mask.any(0)],
        "homdel_partner_genes": keys[seg_of == partner_seg], "cn_only": keys[sorted(cn_only)],
    }
    return ds, truth
