"""Regenerate every briefing figure into docs/briefing/figures/ (SVG + 150-dpi PNG) plus figures.json.

Usage:  .venv/bin/python scripts/figures/make_figures.py

Read-only use of committed results. Sources:
  TEST  : results/final/summary.tsv, results/final/<run>/<line>_TEST*.json
  TRAIN : results/tuning/*.json, results/m4/*.txt, results/tuning_log.tsv,
          data/processed/bench/*TRAIN*.parquet, data/processed/bench/universe_regulators_*.tsv (TRAIN rows only)
  other : results/split/arm_units.tsv, data/external/annotation/cytoBand.hg38.txt.gz,
          docs/DATA_PROFILE.md §10, git log (commit timestamps)
Never reads the sealed directory or any *_TEST* label file.

Style follows the dataviz method: neutral ink for text, one axis per panel, thin marks, hairline
recessive grid, direct labels, fixed categorical order (slot 1 blue, slot 2 orange, slot 3 aqua;
validated light + dark, adjacent and all-pairs). SVGs carry <title> tooltips on data marks.
"""
from __future__ import annotations

import gzip
import io
import json
import math
import re
import subprocess
from datetime import datetime
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "briefing" / "figures"
FINAL = ROOT / "results" / "final"
TUNING = ROOT / "results" / "tuning"

# ---------------------------------------------------------------- tokens (palette.md, light)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
WASH = "#f0efec"          # diverging-midpoint gray, used as a neutral panel wash
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"   # categorical slots 1-3 (fixed order)
BLUE_100, BLUE_150 = "#cde2fb", "#b7d3f6"            # sequential blue steps (fills only)
FAINT = "#b9b7ae"         # de-emphasis gray for rejected / context marks
CATEGORICAL = [BLUE, ORANGE, AQUA]

plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"],
    "font.size": 9,
    "text.color": INK,
    "axes.labelcolor": INK2,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlecolor": INK, "axes.titlepad": 10,
    "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round",
    "svg.fonttype": "none",
    "svg.hashsalt": "depmaptr-figures",
})

# ---------------------------------------------------------------- tooltip + save helpers
_TIPS: dict[str, str] = {}


def tip(artist, text: str):
    """Attach an SVG <title> tooltip to an artist (injected after save)."""
    gid = f"tt{len(_TIPS):04d}"
    artist.set_gid(gid)
    _TIPS[gid] = text
    return artist


def _inject_titles(svg: str) -> str:
    def rep(m):
        gid = m.group(1)
        return f'{m.group(0)}<title>{escape(_TIPS[gid])}</title>'
    return re.sub(r'<g id="(tt\d{4})"[^>]*>', rep, svg)


MANIFEST: list[dict] = []


def save(fig, name: str, title: str, caption: str, sources: list[str]):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png", dpi=150)
    buf = io.StringIO()
    fig.savefig(buf, format="svg", metadata={"Date": None})
    (OUT / f"{name}.svg").write_text(_inject_titles(buf.getvalue()))
    plt.close(fig)
    _TIPS.clear()
    MANIFEST.append({"id": name, "svg": f"{name}.svg", "png": f"{name}.png", "title": title,
                     "caption": caption, "sources": sources})
    print(f"wrote {name}")


def style_axis(ax, grid_axis="x"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.8, linestyle="-", zorder=0)
    ax.set_axisbelow(True)


def log_or_axis(ax, ticks, axis="x"):
    fmt = FuncFormatter(lambda v, _: f"{v:g}")
    a = ax.xaxis if axis == "x" else ax.yaxis
    (ax.set_xscale if axis == "x" else ax.set_yscale)("log")
    a.set_major_locator(FixedLocator(ticks))
    a.set_major_formatter(fmt)
    a.set_minor_locator(NullLocator())


def dot(ax, x, y, color, size=6.5, hollow=False, z=4, marker="o", **kw):
    """Marker with a 2px surface ring (marks-and-anatomy)."""
    (h,) = ax.plot([x], [y], marker=marker, ms=size, mfc=SURFACE if hollow else color,
                   mec=color if hollow else SURFACE, mew=1.6 if hollow else 1.4, ls="none", zorder=z, **kw)
    return h


# ---------------------------------------------------------------- data access
def final_json(run: str, line: str = "hct116", subset: str | None = None) -> dict:
    name = f"{line}_TEST" + (f"_{subset}" if subset else "") + ".json"
    return json.loads((FINAL / run / name).read_text())["metrics"]


def pooled_summary() -> pd.DataFrame:
    d = pd.read_csv(FINAL / "summary.tsv", sep="\t")
    return d[(d.metric_set == "pooled") & d.subset.isna()].set_index(["run", "line"])


def parse_objs(path: Path) -> list[dict]:
    """All JSON objects in a concatenated-JSON file (same parsing as scripts/summarise_compare.py)."""
    txt = path.read_text()
    dec, i, objs = json.JSONDecoder(), 0, []
    while (j := txt.find("{", i)) >= 0:
        try:
            o, e = dec.raw_decode(txt[j:])
            objs.append(o)
            i = j + e
        except json.JSONDecodeError:
            i = j + 1
    return objs


def parse_compare(path: Path) -> list[dict]:
    return [o for o in parse_objs(path) if "candidate" in o]


def train_or(run_id: str, line: str = "hct116") -> float:
    d = pd.read_csv(ROOT / "results" / "tuning_log.tsv", sep="\t")
    d = d[(d.kind == "evaluate") & (d.run_id == run_id) & (d.line == line) & (d.split == "TRAIN")]
    return float(json.loads(d.iloc[-1].metrics)["primary_or"])


# ================================================================= F1 pipeline
def fig_pipeline():
    fig = plt.figure(figsize=(11, 5.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 110)
    ax.set_ylim(0, 54)
    ax.axis("off")

    def box(x, y, w, h, fill=SURFACE, edge=AXIS, lw=1.0):
        p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=fill, ec=edge, lw=lw,
                           zorder=2)
        ax.add_patch(p)
        return p

    def badge(x, y, label, color=BLUE):
        ax.add_patch(plt.Circle((x, y), 1.5, fc=color, ec=SURFACE, lw=1.5, zorder=4))
        ax.text(x, y - 0.05, label, ha="center", va="center", color="#ffffff", fontsize=9, fontweight="bold",
                zorder=5)

    def arrow(x0, y0, x1, y1, color=MUTED):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>,head_length=5,head_width=3",
                                     color=color, lw=1.2, zorder=1, shrinkA=0, shrinkB=0))

    top, h = 31.5, 20.5
    cols = [(1.5, 18.0), (23.0, 20.0), (46.5, 19.5), (69.5, 19.5), (92.5, 16.0)]
    x, w = cols[0]
    box(x, top, w, h, fill=WASH, edge=WASH)
    ax.text(x + 1.4, top + h - 2.0, "Inputs", fontweight="bold", va="top", fontsize=10)
    ax.text(x + 1.4, top + h - 5.2, "DepMap 26Q3\n1,436 cell lines", va="top", color=INK, fontsize=8.6,
            linespacing=1.35)
    ax.text(x + 1.4, top + h - 10.2,
            "• gene copy number (WGS)\n• mRNA, log2(TPM+1)\n• covariates: lineage,\n  ploidy/WGD, MSI, TP53,\n"
            "  sex, media, RNA batch", va="top", color=INK2, fontsize=8, linespacing=1.35)

    steps = [
        ("1", "Trans association", "for every pair A→B on\ndifferent chromosomes",
         "B's mRNA ~ A's CN, adjusted\nfor B's own CN, arm CN of all\nother chromosomes (A's own\n"
         "chromosome left out) and\ncovariates → t-statistic"),
        ("2", "Cis-instrument gate", "keep regulator A only if\nits own mRNA follows its CN",
         "first-stage F ≥ 100\n(positive cis t ≥ 10); other\nregulators get weight 0"),
        ("3", "Co-segment credit", "neighbouring genes share\nthe same CN; share credit",
         "SuSiE L=1 posterior within\neach co-segment (r ≥ 0.95),\nprior ∝ exp(t_cis / 2)"),
    ]
    for (x, w), (n, ttl, what, how) in zip(cols[1:4], steps):
        p = box(x, top, w, h)
        tip(p, f"Step {n}: {ttl}")
        badge(x + 2.6, top + h - 2.8, n)
        ax.text(x + 4.8, top + h - 2.8, ttl, fontweight="bold", va="center", fontsize=9.6)
        ax.text(x + 1.4, top + h - 5.8, what, va="top", color=INK, fontsize=8.4, linespacing=1.35)
        ax.text(x + 1.4, top + h - 10.6, how, va="top", color=INK2, fontsize=7.9, linespacing=1.35)

    x, w = cols[4]
    p = box(x, top, w, h, fill=BLUE_100, edge=BLUE_150)
    tip(p, "Output of the selected CN-only tool (M1 + G1 + A3)")
    ax.text(x + 1.4, top + h - 2.0, "Output", fontweight="bold", va="top", fontsize=10)
    ax.text(x + 1.4, top + h - 5.2, "ranked, signed trans\npredictions", va="top", fontsize=8.6,
            linespacing=1.35)
    ax.text(x + 1.4, top + h - 10.2,
            "score = t × gate × credit\nsign → predicted direction\nafter knockdown\n≈2.2×10⁸ pairs in ~67 s\n"
            "(laptop CPU)", va="top", color=INK2, fontsize=7.9, linespacing=1.35)

    ymid = top + h / 2
    for (x0, w0), (x1, _) in zip(cols[:-1], cols[1:]):
        arrow(x0 + w0 + 0.4, ymid, x1 - 0.5, ymid)

    # optional context layer, beneath the output
    cx, cw, cy, ch = 96.5, 12.0, 15.5, 12.5
    p = box(cx, cy, cw, ch, fill=SURFACE, edge=ORANGE)
    tip(p, "Context layer: reported variant, not the selected method (its gain survives the pair-shuffle control)")
    ax.text(cx + 1.2, cy + ch - 1.6, "Optional\ncontext layer", fontweight="bold", va="top", fontsize=8.8,
            linespacing=1.2)
    ax.text(cx + 1.2, cy + ch - 6.2, "re-rank by how\ndetectable target B\nis in the context line", va="top",
            color=INK2, fontsize=7.6, linespacing=1.3)
    arrow(cx + cw / 2, top - 0.3, cx + cw / 2, cy + ch + 0.4)

    # evaluation band
    ex, ew, ey, eh = 1.5, 107.0, 1.5, 11.5
    box(ex, ey, ew, eh, fill=WASH, edge=WASH)
    ax.text(ex + 1.4, ey + eh - 1.6, "Evaluation against CRISPRi Perturb-seq", fontweight="bold", va="top",
            fontsize=9.6)
    ax.text(ex + 1.4, ey + eh - 4.6,
            "Arc Institute PIE X-Atlas/Orion DE tables; HCT116 primary, HEK293T secondary\n"
            "positive pair = knocking down A changes B's mRNA (BH p_adj < 0.05)\n"
            "held out: regulators on 14 whole chromosome arms, labels sealed before tuning",
            va="top", color=INK2, fontsize=8, linespacing=1.4)
    ax.text(ex + 50, ey + eh - 4.6,
            "odds ratio of hits among the top 0.1% of pairs\n(regulator-label permutation p, arm-jackknife CI)\n"
            "sign concordance: predicted vs observed direction",
            va="top", color=INK2, fontsize=8, linespacing=1.4)
    arrow(cx + cw / 2, cy - 0.3, cx + cw / 2, ey + eh + 0.4)
    arrow(93.5, top - 0.3, 93.5, ey + eh + 0.4)
    save(fig, "pipeline", "The final method",
         "The selected CN-only tool: covariate-adjusted trans association with leave-own-chromosome-out arm "
         "adjustment, a cis-instrument gate and co-segment attribution, with an optional context layer; both "
         "outputs are evaluated against held-out CRISPRi Perturb-seq.",
         ["REPORT.md §1–2", "docs/MODEL.md", "results/final/selected.yaml"])


# ================================================================= F2 benchmark structure
def _parse_profile_row(text: str, label_regex: str) -> list[str]:
    m = re.search(rf"^\|\s*{label_regex}\s*\|(.*)\|\s*$", text, flags=re.M)
    if not m:
        raise SystemExit(f"DATA_PROFILE row not found: {label_regex}")
    return [c.strip().strip("*") for c in m.group(1).split("|")]


def fig_benchmark_structure():
    lines = [("hct116", "HCT116", BLUE), ("hek293t", "HEK293T", ORANGE)]
    curves = {}
    for key, _, _ in lines:
        u = pd.read_csv(ROOT / f"data/processed/bench/universe_regulators_{key}.tsv", sep="\t", usecols=["regulator", "split"])
        u = u[u.split == "TRAIN"]
        p = pd.read_parquet(ROOT / f"data/processed/bench/positives_{key}_TRAIN.parquet",
                            columns=["regulator", "pos_primary"])
        p = p[p.pos_primary]
        c = p.groupby("regulator").size().reindex(u.regulator, fill_value=0).sort_values(ascending=False).values
        curves[key] = c

    prof = (ROOT / "docs" / "DATA_PROFILE.md").read_text()
    mean_row = _parse_profile_row(prof, r"mean trans hits \(0\.05\): own gene expressed vs not")
    share_row = _parse_profile_row(prof, r"share of all trans hits from non-expressed targets")
    expr_row = _parse_profile_row(prof, r"targets whose own gene is tested \(expressed\)")
    tot_row = _parse_profile_row(prof, r"targets / features")
    means = [tuple(float(v) for v in s.split(" vs ")) for s in mean_row]
    shares = [s for s in share_row]
    n_expr = [s.split(" ")[0] for s in expr_row]
    n_tot = [s.split(" / ")[0] for s in tot_row]

    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.3), gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.32})
    fig.subplots_adjust(left=0.065, right=0.975, top=0.86, bottom=0.17)
    # (a) concentration curve (log x so the steep head is visible)
    xd = np.logspace(math.log10(0.02), 2, 200)
    a.plot(xd, xd, color=AXIS, lw=1, zorder=1)
    a.text(4.2, 2.2, "equal share", color=MUTED, fontsize=8, rotation=0)
    for key, lab, col in lines:
        c = curves[key]
        n = len(c)
        xs = np.arange(1, n + 1) / n * 100
        ys = np.cumsum(c) / c.sum() * 100
        h, = a.plot(xs, ys, color=col, lw=2, zorder=3)
        k10 = int(math.ceil(0.10 * n))
        k1 = int(math.ceil(0.01 * n))
        tip(h, f"{lab} TRAIN: {n:,} perturbations, {c.sum():,} trans hits; top 1% → {ys[k1 - 1]:.0f}%, "
               f"top 10% → {ys[k10 - 1]:.0f}%; {np.mean(c == 0) * 100:.0f}% have no hit")
        curves[key] = (n, c.sum(), ys[k1 - 1], ys[k10 - 1], np.mean(c == 0) * 100)
    n1, h1, s1_1, s1_10, z1 = curves["hct116"]
    n2, h2, s2_1, s2_10, z2 = curves["hek293t"]
    for xv, v1, v2 in ((1, s1_1, s2_1), (10, s1_10, s2_10)):
        a.axvline(xv, color=GRID, lw=0.8, zorder=0)
        tip(dot(a, xv, v1, BLUE, size=6.5), f"HCT116: top {xv}% of perturbations → {v1:.0f}% of hits")
        tip(dot(a, xv, v2, ORANGE, size=6.5), f"HEK293T: top {xv}% of perturbations → {v2:.0f}% of hits")
    a.text(1.12, 30, f"top 1% of perturbations:\n{s1_1:.0f}% (HCT116), {s2_1:.0f}% (HEK293T)\nof all trans hits",
           fontsize=8, color=INK2, va="top")
    a.text(11.2, 80, f"top 10%:\n{s1_10:.0f}% / {s2_10:.0f}%", fontsize=8, color=INK2, va="top")
    for (lab, col, nn, hh), yl in zip((("HCT116", BLUE, n1, h1), ("HEK293T", ORANGE, n2, h2)), (95, 80)):
        a.plot([0.03, 0.06], [yl, yl], color=col, lw=2)
        a.text(0.07, yl, lab, fontsize=8.3, color=INK, va="center")
        a.text(0.07, yl - 6, f"{nn:,} perturbations, {hh:,} hits", fontsize=7.6, color=INK2, va="center")
    log_or_axis(a, [0.1, 1, 10, 100])
    a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}%"))
    a.set_xlim(0.02, 100)
    a.set_ylim(0, 102)
    a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}%"))
    a.set_xlabel("top share of perturbations, ranked by number of trans hits (log scale)")
    a.set_ylabel("share of all trans hits captured")
    a.set_title("a  A few perturbations carry almost all trans hits")
    style_axis(a, "y")
    fig.text(0.065, 0.025, f"(a) TRAIN regulators in the evaluation universe; {z1:.0f}% (HCT116) and {z2:.0f}% "
             "(HEK293T) have no hit at all.  (b) All perturbations in the full DE tables (docs/DATA_PROFILE.md §10).",
             fontsize=7.5, color=MUTED)

    # (b) dumbbell: mean hits per perturbation, own gene expressed vs not
    ys = [1, 0]
    for (key, lab, col), y, (m_e, m_n), sh, ne, nt in zip(lines, ys, means, shares, n_expr, n_tot):
        b.plot([m_n, m_e], [y, y], color=AXIS, lw=2, zorder=2, solid_capstyle="butt")
        tip(dot(b, m_e, y, col, size=8.5), f"{lab}: perturbed gene expressed ({ne} of {nt}): mean {m_e} trans hits")
        tip(dot(b, m_n, y, col, size=8.5, hollow=True),
            f"{lab}: perturbed gene not expressed: mean {m_n} trans hits; {sh} of all trans hits")
        b.text(m_e, y + 0.17, f"{m_e:g}", ha="center", fontsize=8.5, color=INK)
        b.text(m_n, y + 0.17, f"{m_n:g}", ha="center", fontsize=8.5, color=INK)
        b.text(62, y - 0.02, f"{sh} of all trans hits come\nfrom non-expressed perturbations",
               fontsize=8, color=INK2, va="center")
    b.set_yticks(ys, [f"HCT116\n{n_expr[0]} of {n_tot[0]} expressed", f"HEK293T\n{n_expr[1]} of {n_tot[1]} expressed"])
    b.tick_params(axis="y", length=0)
    b.set_ylim(-0.6, 1.75)
    b.set_xlim(0, 100)
    b.set_xticks([0, 10, 20, 30, 40, 50, 60])
    b.spines["bottom"].set_bounds(0, 60)
    b.set_xlabel("mean significant trans hits per perturbation")
    b.set_title("b  Non-expressed perturbed genes still produce hits")
    style_axis(b, None)
    b.spines["left"].set_visible(False)
    for xg in [10, 20, 30, 40, 50, 60]:
        b.axvline(xg, color=GRID, lw=0.8, zorder=0, ymax=0.83)
    # legend (two marker kinds)
    dot(b, 2, 1.58, INK2, size=7)
    b.text(4.6, 1.58, "perturbed gene expressed in the line", va="center", fontsize=8, color=INK2)
    dot(b, 2, 1.40, INK2, size=7, hollow=True)
    b.text(4.6, 1.40, "not expressed (knockdown cannot act; expected ≈ null)", va="center", fontsize=8, color=INK2)
    save(fig, "benchmark_structure", "Why the benchmark is hard",
         f"Perturb-seq trans hits are concentrated in a few perturbations (top 10% carry {s1_10:.0f}% of HCT116 "
         f"and {s2_10:.0f}% of HEK293T TRAIN hits), and in HEK293T most hits come from perturbations whose gene is "
         "not even expressed, so HCT116 is the primary line.",
         ["data/processed/bench/positives_hct116_TRAIN.parquet", "data/processed/bench/positives_hek293t_TRAIN.parquet",
          "data/processed/bench/universe_regulators_hct116.tsv (TRAIN rows)",
          "data/processed/bench/universe_regulators_hek293t.tsv (TRAIN rows)", "docs/DATA_PROFILE.md §10–11"])


# ================================================================= F3 split ideogram
def fig_split():
    arms = pd.read_csv(ROOT / "results/split/arm_units.tsv", sep="\t").set_index("arm")
    cyto = pd.read_csv(gzip.open(ROOT / "data/external/annotation/cytoBand.hg38.txt.gz", "rt"), sep="\t",
                       header=None, names=["chrom", "start", "end", "band", "stain"])
    chroms = [str(i) for i in range(1, 23)] + ["X"]
    geo = {}
    for c in chroms:
        sub = cyto[cyto.chrom == f"chr{c}"]
        cen = sub[(sub.stain == "acen") & sub.band.str.startswith("p")].end.max()
        geo[c] = (sub.end.max() / 1e6, cen / 1e6)
    fig, ax = plt.subplots(figsize=(11, 4.6))
    fig.subplots_adjust(left=0.03, right=0.99, top=0.80, bottom=0.08)
    maxlen = max(v[0] for v in geo.values())
    w = 0.42
    col = {"TRAIN": BLUE, "TEST": ORANGE}
    for i, c in enumerate(chroms):
        length, cen = geo[c]
        for arm, y0, y1 in ((f"{c}p", 0, cen), (f"{c}q", cen, length)):
            gap = 0.9  # Mb-equivalent surface gap at the centromere
            lo, hi = (y0, y1 - gap / 2) if arm.endswith("p") else (y0 + gap / 2, y1)
            if arm in arms.index:
                r = arms.loc[arm]
                fc = col[r.split]
                txt = f"{arm}: {r.split} · {r.n_hct116:,} HCT116 / {r.n_hek293t:,} HEK293T regulators"
            else:
                fc = WASH
                txt = f"{arm}: no regulators in the universe (acrocentric p arm)"
            rect = Rectangle((i - w / 2, -hi), w, hi - lo, fc=fc, ec="none", zorder=2)
            ax.add_patch(rect)
            tip(rect, txt)
            if arm in arms.index:
                n = int(arms.loc[arm].n_hct116)
                if arm.endswith("p"):
                    ax.text(i, -lo + 4, f"{n}", ha="center", va="bottom", fontsize=7, color=INK2)
                else:
                    ax.text(i, -hi - 4, f"{n}", ha="center", va="top", fontsize=7, color=INK2)
        ax.text(i, 34, c, ha="center", va="bottom", fontsize=9, color=INK, fontweight="bold")
    ax.set_xlim(-0.7, len(chroms) - 0.3)
    ax.set_ylim(-maxlen - 22, 44)
    ax.axis("off")
    tr = arms[arms.split == "TRAIN"]
    te = arms[arms.split == "TEST"]
    fig.text(0.03, 0.94, "Held-out design: whole chromosome arms", fontsize=10.5, fontweight="bold", color=INK)
    fig.text(0.03, 0.885, "numbers = regulators (perturbed genes) per arm in the HCT116 evaluation universe; "
             "p arm above the gap, q arm below; gray = acrocentric p arms with no regulators", fontsize=8, color=INK2)
    lx, ly = 0.60, 0.20
    for k, (c, txt) in enumerate([(BLUE, f"TRAIN: {len(tr)} arms, {tr.n_hct116.sum():,} regulators (all tuning)"),
                                  (ORANGE, f"TEST: {len(te)} arms, {te.n_hct116.sum():,} regulators (sealed, used once)")]):
        yy = ly - k * 0.06
        fig.patches.append(Rectangle((lx, yy - 0.015), 0.014, 0.03, fc=c, ec="none", transform=fig.transFigure))
        fig.text(lx + 0.02, yy, txt, va="center", fontsize=8.8, color=INK)
    save(fig, "split", "Held-out design",
         f"Regulators were split by whole chromosome arm: {len(tr)} TRAIN arms ({tr.n_hct116.sum():,} HCT116 "
         f"regulators) for all tuning and {len(te)} TEST arms ({te.n_hct116.sum():,} regulators) whose Perturb-seq "
         "labels were sealed until the single final run.",
         ["results/split/arm_units.tsv", "data/external/annotation/cytoBand.hg38.txt.gz"])


# ================================================================= F4 tuning path
TIER_FILES = [
    ("1 arm handling", ROOT / "results/m4/tier1_compare.txt"),
    ("2 G1 cis gate", TUNING / "tier2_G1_compare.json"),
    ("3 G2 context expr.", TUNING / "tier3_G2_compare.json"),
    ("4 A3 attribution", TUNING / "tier4_A3_compare.json"),
    ("5–7 transforms", [TUNING / "tier5_tstd_compare.json", TUNING / "tier6_hinge_compare.json",
                        TUNING / "tier7_transforms_compare.json"]),
    ("D6 post-hoc", [TUNING / "d6_1_I5_compare.json", TUNING / "d6_2_I4_compare.json",
                     TUNING / "d6_3_I3_I6_compare.json"]),
    ("D6 context layer", TUNING / "d6_4_I1_compare.json"),
]
SHORT = {
    "focal": "focal arms", "loao": "leave-one-arm-out", "cnpcs20": "CN PCs",
    "g1_t10_w0": "τ=10", "g1_t30_w0": "τ=30", "g1_t100_w0": "τ=100",
    "g1_t100_w0__g2_thr1_w0": "G2 (4 settings)",
    "g1_t100_w0__a3_r095_s2_w0": "ρ .95, s 2", "g1_t100_w0__a3_r098_s2_w0": "ρ .98",
    "g1_t100_w0__a3_r095_s5_w0": "s 5", "g1_t100_w0__a3_r095_s2_w025": "w_co .25",
    "g1_t100_w0__a3_r095_s2_w05": "w_co .5",
    "g1_t100_w0__a3_r095_s2_w0__tstd": "target std.", "g1_t100_w0__a3_r095_s2_w0__hinge": "hinge",
    "g1_t100_w0__a3_r095_s2_w0__floor5": "CN floor", "g1_t100_w0__a3_r095_s2_w0__rint": "RINT",
    "g1_t100_w0__a3_r095_s2_w0__i5_comp": "compositional", "g1_t100_w0__i4_r090": "ρ .90",
    "g1_t100_w0__a3_r095_s2_w0__i3_ivshrink": "IV-shrink", "g1_t100_w0__a3_r095_s2_w0__i6_linrep": "lineage rep.",
    "g1_t100_w0__a3_r095_s2_w0__i1_fuse_g1": "γ=1", "g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2": "γ=2",
}
ADOPTED = {"g1_t100_w0", "g1_t100_w0__a3_r095_s2_w0"}
CONTEXT = {"g1_t100_w0__a3_r095_s2_w0__i1_fuse_g1", "g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2"}


def fig_tuning_path():
    or_m1, or_g1, or_sel = train_or("M1"), train_or("g1_t100_w0"), train_or("g1_t100_w0__a3_r095_s2_w0")
    c1 = [o for o in parse_objs(TUNING / "d6_I1_c1.json")  # exact-mode C1 for the γ=2 context layer
          if o.get("run", "").endswith("i1_fuse_g2") and o.get("c1_mode", "").startswith("pair_score")][0]
    or_ctx, or_ctx_c1 = c1["or_run"], c1["or_c1"]
    refs = {b: train_or(b) for b in ["B5", "B4e", "B2", "B0"]}

    tiers = []
    for name, f in TIER_FILES:
        objs = []
        for p in (f if isinstance(f, list) else [f]):
            objs += parse_compare(p)
        seen, rows = set(), []
        for o in objs:
            cand = o["candidate"]
            if cand.endswith("_w05") and "a3" not in cand and "g2" not in cand:
                continue  # G1 w=0.5 is numerically identical to w=0 (see tier-2 file)
            if "__g2_" in cand:
                if "g2" in seen:
                    continue
                seen.add("g2")
                cand = "g1_t100_w0__g2_thr1_w0"
            rows.append(dict(cand=cand, d=o["delta_log_or"], f1=o["delta_log_or_folds"]["fold1"],
                             f2=o["delta_log_or_folds"]["fold2"], lb=o["boot_delta_log_or_lb90"],
                             rule=o["adopt_rule_8_2"], inc=o["incumbent"]))
        tiers.append((name, rows))

    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 6.6), gridspec_kw={"width_ratios": [1, 1.0], "wspace": 0.62})
    fig.subplots_adjust(left=0.065, right=0.975, top=0.90, bottom=0.135)

    # (a) adopted path, log-scale OR
    xs = [0, 1, 2]
    ys = [or_m1, or_g1, or_sel]
    for name, v in refs.items():
        a.axhline(v, color=GRID, lw=1, zorder=1)
        lab = {"B5": "B5 hub prior", "B4e": "B4e target expr.", "B2": "B2 co-expression", "B0": "B0 random"}[name]
        if name == "B0":
            a.text(3.65, v * 1.04, f"{lab} {v:.2f}", fontsize=7.5, color=MUTED, va="bottom", ha="right")
        else:
            a.text(-0.3, v * (0.965 if name == "B2" else 1.035), f"{lab} {v:.1f}", fontsize=7.5, color=MUTED,
                   va="top" if name == "B2" else "bottom", ha="left")
    a.plot(xs, ys, color=BLUE, lw=2, zorder=3)
    for x, y, t in zip(xs, ys, ["M1 core model", "+ G1 adopted", "+ A3 adopted"]):
        tip(dot(a, x, y, BLUE, size=8), f"{t}: TRAIN OR@0.1% = {y:.2f}")
        a.text(x, y * 1.22, f"{y:.2f}", ha="center", fontsize=9, color=INK, fontweight="bold")
    a.plot([2, 3], [or_sel, or_ctx], color=FAINT, lw=1.2, zorder=2)
    tip(dot(a, 3, or_ctx, ORANGE, size=8),
        f"context layer (γ=2): TRAIN OR {or_ctx:.1f}; pair-shuffle control C1 {or_ctx_c1:.1f}")
    a.text(3, or_ctx * 1.22, f"{or_ctx:.1f}", ha="center", fontsize=9, color=INK, fontweight="bold")
    a.text(2.93, or_ctx * 0.56, f"C1 shuffle {or_ctx_c1:.1f}\n→ not selected", ha="right", fontsize=7.5,
           color=INK2, va="top")
    log_or_axis(a, [1, 2, 5, 10, 20, 50], "y")
    a.set_ylim(0.8, 60)
    a.set_xlim(-0.35, 3.7)
    a.set_xticks([0, 1, 2, 3], ["M1", "+ G1\ncis gate", "+ A3\nattribution\n(selected)", "context layer\n(post-hoc)"])
    a.tick_params(axis="x", length=0)
    a.set_ylabel("TRAIN odds ratio, top 0.1% (log scale)")
    a.set_title("a  Adopted path (HCT116 TRAIN)")
    style_axis(a, "y")

    # (b) every candidate, one row each, Δ log-OR vs the incumbent at that step
    y = 0
    for name, rows in tiers:
        b.text(-0.03, y, name, transform=b.get_yaxis_transform(), ha="right", va="center", fontsize=8.6,
               fontweight="bold", color=INK)
        y += 1
        for q in sorted(rows, key=lambda q: -q["d"]):
            if q["cand"] in ADOPTED:
                col = BLUE
            elif q["cand"] in CONTEXT:
                col = ORANGE
            else:
                col = FAINT
            b.plot([0, q["d"]], [y, y], color=GRID, lw=2, zorder=1, solid_capstyle="butt")
            h = dot(b, q["d"], y, col, size=7 if col != FAINT else 6)
            status = ("adopted" if q["cand"] in ADOPTED else
                      "passes rule; reported, not selected (fails C1)" if q["cand"] in CONTEXT else
                      "passes rule but a larger gain was adopted" if q["rule"] else "rejected")
            lab = SHORT.get(q["cand"], q["cand"])
            tip(h, f"{lab}: Δlog-OR {q['d']:+.2f} vs {q['inc']} "
                   f"(folds {q['f1']:+.2f}/{q['f2']:+.2f}; bootstrap lb90 {q['lb']:+.2f}) — {status}")
            b.text(-0.03, y, lab, transform=b.get_yaxis_transform(), ha="right", va="center", fontsize=8,
                   color=INK if col != FAINT else INK2)
            xt = max(q["d"], 0) + 0.04
            b.text(xt, y, f"{q['d']:+.2f}" + ("  adopted" if col == BLUE else "  reported" if q["cand"].endswith("g2") else ""),
                   va="center", fontsize=7.6, color=INK if col != FAINT else INK2)
            y += 1
        y += 0.4
    b.axvline(0, color=INK2, lw=1, zorder=1)
    b.set_yticks([])
    b.set_ylim(y - 0.5, -0.8)
    b.set_xlim(-0.3, 2.2)
    b.set_xticks([0, 0.5, 1.0, 1.5])
    b.set_xlabel("Δ log odds ratio vs the incumbent at that step")
    b.set_title("b  Every candidate tried, in order")
    style_axis(b, "x")
    b.spines["left"].set_visible(False)
    fig.text(0.065, 0.015, "Adoption rule (pre-registered): Δ log-OR > 0 in both TRAIN arm-folds and paired-bootstrap "
             "90% lower bound > 0. Tiers 1–2 compare against M1, tiers 3–4 against M1 + G1, later tiers against M1 + G1 + A3.\n"
             "G1 τ=10/30 and A3 w_co .25 also passed the rule; the largest gain was adopted. G1 weight w=0.5 gave identical results.",
             fontsize=7.5, color=MUTED)
    save(fig, "tuning_path", "How the method was chosen",
         f"On HCT116 TRAIN, the cis-instrument gate (G1) and co-segment attribution (A3) raised the top-0.1% odds "
         f"ratio from {or_m1:.2f} to {or_sel:.2f}; every other candidate was rejected, and the context layer "
         f"({or_ctx:.1f}) was reported but not selected because a within-regulator pair shuffle kept its gain.",
         ["results/tuning_log.tsv", "results/m4/tier1_compare.txt", "results/tuning/tier2_G1_compare.json",
          "results/tuning/tier3_G2_compare.json", "results/tuning/tier4_A3_compare.json",
          "results/tuning/tier5_tstd_compare.json", "results/tuning/tier6_hinge_compare.json",
          "results/tuning/tier7_transforms_compare.json", "results/tuning/d6_1_I5_compare.json",
          "results/tuning/d6_2_I4_compare.json", "results/tuning/d6_3_I3_I6_compare.json",
          "results/tuning/d6_4_I1_compare.json", "results/tuning/d6_I1_c1.json"])


# ================================================================= F5 TEST forest
FOREST = [
    ("Method", [("g1_t100_w0__a3_r095_s2_w0", "Selected CN-only tool (M1 + G1 + A3)")]),
    ("Ablations / variants", [("a3_r095_s2_w0", "minus cis gate G1"),
                              ("focal__g1_t100_w0__a3_r095_s2_w0", "focal arm handling"),
                              ("g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2", "+ context layer")]),
    ("Baselines", [("B0", "B0 random"), ("B1", "B1 naive CN–mRNA correlation"),
                   ("B2", "B2 mRNA co-expression"), ("B4", "B4 target expression SD"),
                   ("B4e", "B4e target expression in HCT116"), ("B5", "B5 DepMap hub prior"),
                   ("B6", "B6 arm-level CN only")]),
    ("Controls", [("NC1", "NC1 permuted CN"), ("C2_g2", "C2 context layer on permuted CN")]),
]


def fmt_or(v: float) -> str:
    return f"{v:.2f}" if v < 10 else f"{v:.1f}" if v < 100 else f"{v:.0f}"


def fig_test_results():
    s = pooled_summary()
    rows, y, ypos, headers = [], 0, [], []
    for g, items in FOREST:
        headers.append((g, y))
        y += 1
        for run, lab in items:
            rows.append((run, lab, y))
            y += 1
        y += 0.35
    ymax = y
    fig = plt.figure(figsize=(11, 6.4))
    a = fig.add_axes([0.27, 0.13, 0.46, 0.75])
    t = fig.add_axes([0.737, 0.13, 0.088, 0.75])
    b = fig.add_axes([0.845, 0.13, 0.14, 0.75])
    for ax in (a, t, b):
        ax.set_ylim(ymax - 0.2, -0.8)
    t.axis("off")
    log_or_axis(a, [0.3, 1, 3, 10, 30, 100])
    a.set_xlim(0.25, 250)
    log_or_axis(b, [0.01, 0.1, 1, 10, 100])
    b.set_xlim(0.008, 500)
    for ax in (a, b):
        style_axis(ax, "x")
        ax.spines["left"].set_visible(False)
        ax.set_yticks([])
        ax.axvline(1, color=INK2, lw=1, zorder=1)
    for g, yy in headers:
        fig.text(0.015, a.transData.transform((1, yy))[1] / fig.bbox.height, g, fontsize=9, fontweight="bold",
                 va="center", color=INK)
    for run, lab, yy in rows:
        sel = run == "g1_t100_w0__a3_r095_s2_w0"
        col = BLUE if sel else INK2
        yf = a.transData.transform((1, yy))[1] / fig.bbox.height
        fig.text(0.03, yf, lab, fontsize=8.6, va="center", color=INK, fontweight="bold" if sel else "normal")
        for ax, line in ((a, "hct116"), (b, "hek293t")):
            r = s.loc[(run, line)]
            o, lo, hi, pp = r["OR@0.1%"], r["jk_CI95_lo"], r["jk_CI95_hi"], r["perm_p"]
            ln, = ax.plot([lo, hi], [yy, yy], color=col, lw=1.6 if sel else 1.2, zorder=2,
                          solid_capstyle="butt")
            m = dot(ax, o, yy, col, size=7.5 if sel else 6)
            tip(m, f"{lab}, {line.upper()} TEST: OR {o:.2f} (jackknife 95% CI {lo:.2f}–{hi:.2f}); "
                   f"permutation p = {pp:.3f}")
            tip(ln, f"{lab}: jackknife 95% CI {lo:.2f}–{hi:.2f}")
            if line == "hct116":
                t.text(0.0, yy, fmt_or(o), ha="left", va="center", fontsize=8.3, color=INK,
                       fontweight="bold" if sel else "normal", transform=t.get_yaxis_transform())
                t.text(1.0, yy, f"{fmt_or(lo)}–{fmt_or(hi)}", ha="right", va="center", fontsize=7.6, color=INK2,
                       transform=t.get_yaxis_transform())
            if sel and line == "hct116":
                ax.annotate(f"pre-registered permutation p = {pp:.3f}", xy=(hi, yy), xytext=(hi * 1.35, yy),
                            fontsize=8, color=INK2, va="center")
    t.text(0.0, -0.75, "OR", fontsize=8, color=INK2, transform=t.get_yaxis_transform(), va="center")
    t.text(1.0, -0.75, "95% CI", fontsize=8, color=INK2, ha="right", transform=t.get_yaxis_transform(),
           va="center")
    a.set_title("HCT116 TEST (primary)", pad=8)
    b.set_title("HEK293T TEST", pad=8)
    a.set_xlabel("odds ratio of Perturb-seq hits in top 0.1% of pairs (log scale)")
    b.set_xlabel("odds ratio (log scale)")
    fig.text(0.015, 0.955, "Held-out (TEST) enrichment, with arm-jackknife 95% CIs", fontsize=10.5,
             fontweight="bold", color=INK)
    fig.text(0.015, 0.025, "Reference line at OR = 1 (no enrichment). B3 (cis instrument × CN SD) is omitted: "
             "1 hit in the top 0.1%, CI unbounded. HCT116: 14 arms, 1,817 regulators, 17.7 M pairs.",
             fontsize=7.5, color=MUTED)
    sel = s.loc[("g1_t100_w0__a3_r095_s2_w0", "hct116")]
    save(fig, "test_results", "Held-out TEST results",
         f"On HCT116 TEST the selected CN-only tool reaches OR {sel['OR@0.1%']:.2f} (jackknife CI "
         f"{sel['jk_CI95_lo']:.2f}–{sel['jk_CI95_hi']:.1f}; pre-registered permutation p = {sel['perm_p']:.3f}), "
         "numerically below the target-detectability, co-expression and hub baselines (no paired difference significant), "
         "and it does not replicate in HEK293T.",
         ["results/final/summary.tsv (metric_set = pooled)"])


# ================================================================= F6 sign concordance
SIGN_GROUPS = [  # (group, [(run, subset, label, mark colour)])
    ("Method", [("g1_t100_w0__a3_r095_s2_w0", None, "Selected CN-only tool", BLUE)]),
    ("Sensitivity", [("g1_t100_w0__a3_r095_s2_w0__S1", None, "S1 sample set (999 lines)", INK2),
                     ("g1_t100_w0__a3_r095_s2_w0", "no_split_sibling", "no-split-sibling regulators", INK2)]),
    ("Negative controls (permuted CN)", [("C2_g2", None, "C2: context fusion on permuted CN", ORANGE),
                                         ("NC1", None, "NC1: permuted CN", INK2)]),
]


def fig_sign_concordance():
    fig, ax = plt.subplots(figsize=(10, 4.6))
    fig.subplots_adjust(left=0.27, right=0.80, top=0.76, bottom=0.17)
    y, srcs = 0, []
    for g, rows in SIGN_GROUPS:
        ax.text(-0.36, y, g, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8.8,
                fontweight="bold", color=INK)
        y += 1
        for run, sub, lab, col in rows:
            srcs.append(f"results/final/{run}/hct116_TEST{'_' + sub if sub else ''}.json")
            sg = final_json(run, "hct116", sub)["sign"]
            obs, c0, nm, p, nn = (sg["concordance"] * 100, sg["c0"] * 100, sg["sign_perm_null_mean"] * 100,
                                  sg["sign_perm_p"], sg["n_overlap"])
            hl = col != INK2
            ax.plot([nm, obs], [y, y], color=col if hl else AXIS, lw=2, zorder=2, solid_capstyle="butt")
            tip(dot(ax, c0, y, MUTED, size=6.5, hollow=True, marker="s"),
                f"{lab}: pooled independence expectation c0 = {c0:.1f}%")
            tip(dot(ax, nm, y, MUTED, size=7, marker="D"), f"{lab}: within-regulator permutation null mean {nm:.1f}%")
            tip(dot(ax, obs, y, col, size=8.5 if hl else 7),
                f"{lab}: observed concordance {obs:.1f}% (n = {nn} overlapping pairs; sign permutation p = {p:.3f})")
            ax.text(-0.02, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8.6,
                    color=INK)
            ax.text(1.02, y, f"n = {nn}", transform=ax.get_yaxis_transform(), va="center", fontsize=8, color=INK2)
            ax.text(1.12, y, f"p = {p:.3f}", transform=ax.get_yaxis_transform(), va="center", fontsize=8,
                    color=INK if p < 0.05 else INK2, fontweight="bold" if p < 0.05 else "normal")
            if run == "C2_g2":
                ax.text(obs + 1.0, y, "a null control 'passes'", va="center", fontsize=7.8, color=INK2)
            y += 1
        y += 0.3
    ax.set_ylim(y - 0.6, -0.7)
    ax.set_yticks([])
    ax.set_xlim(30, 80)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.set_xlabel("pairs whose predicted knockdown direction matches Perturb-seq (sign concordance)")
    style_axis(ax, "x")
    ax.spines["left"].set_visible(False)
    ax.text(1.02, -0.55, "overlaps", transform=ax.get_yaxis_transform(), fontsize=7.6, color=MUTED, va="center")
    ax.text(1.12, -0.55, "sign perm.", transform=ax.get_yaxis_transform(), fontsize=7.6, color=MUTED, va="center")
    fig.text(0.02, 0.94, "A permuted-CN control passes the same sign test, so the direction result is not established",
             fontsize=10.5, fontweight="bold", color=INK)
    fig.text(0.02, 0.895, "HCT116 TEST; top-0.1% pairs that are Perturb-seq hits. The within-regulator permutation test "
             "is anti-conservative (D9).", fontsize=8, color=INK2)
    lx = 0.27
    for k, (mk, lab, hol, c) in enumerate([("o", "observed", False, INK2),
                                           ("D", "within-regulator null mean", False, MUTED),
                                           ("s", "pooled independence expectation (c0)", True, MUTED)]):
        xx = lx + [0, 0.11, 0.33][k]
        fig.lines.append(plt.Line2D([xx], [0.835], marker=mk, ms=6.5, mfc=SURFACE if hol else c,
                                    mec=c if hol else SURFACE, mew=1.5, transform=fig.transFigure))
        fig.text(xx + 0.012, 0.835, lab, fontsize=8, color=INK2, va="center")
    fig.text(0.02, 0.025, "Untested explanation: targets respond in correlated modules, which the sign permutation ignores; "
             "it was never calibrated on negative controls (REPORT.md §1, finding 2).",
             fontsize=7.5, color=MUTED)
    sg = final_json("g1_t100_w0__a3_r095_s2_w0")["sign"]
    c2 = final_json("C2_g2")["sign"]
    save(fig, "sign_concordance", "Direction of predicted effects: not established",
         f"The selected tool's sign concordance ({sg['concordance'] * 100:.1f}% vs a within-regulator null of "
         f"{sg['sign_perm_null_mean'] * 100:.1f}%, p = {sg['sign_perm_p']:.3f}, n = {sg['n_overlap']}) looks better "
         f"than chance, but the permuted-CN control C2 passes the same test ({c2['concordance'] * 100:.1f}% vs "
         f"{c2['sign_perm_null_mean'] * 100:.1f}%, p = {c2['sign_perm_p']:.3f}, n = {c2['n_overlap']}), so the test is "
         "miscalibrated and the direction result is not established; it is also not significant in S1 or "
         "no-split-sibling regulators.",
         srcs)


# ================================================================= F7 what the signal is
def fig_what_the_signal_is():
    groups = [
        ("CN-only selected tool", [("g1_t100_w0__a3_r095_s2_w0", "real scores", BLUE),
                                   ("g1_t100_w0__a3_r095_s2_w0_C1", "targets shuffled within each regulator (C1)", AQUA),
                                   ("NC1", "CN permuted (NC1)", ORANGE)]),
        ("Context layer", [("g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2", "real scores", BLUE),
                           ("g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1", "targets shuffled within each regulator (C1)", AQUA),
                           ("C2_g2", "CN permuted, same fusion (C2)", ORANGE)]),
    ]
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1.35, 1], "wspace": 0.55})
    fig.subplots_adjust(left=0.20, right=0.975, top=0.80, bottom=0.19)
    y = 0
    yt, ylab = [], []
    for g, items in groups:
        base = None
        a.text(-0.02, y - 0.75, g, transform=a.get_yaxis_transform(), ha="right", fontsize=9, fontweight="bold",
               color=INK, va="center")
        for run, lab, col in items:
            m = final_json(run)
            o = m["primary_or"]
            lo, hi = (math.exp(v) for v in m["jackknife_log_or"]["ci95"])
            base = base or o
            a.plot([lo, hi], [y, y], color=col, lw=1.3, zorder=2, solid_capstyle="butt")
            tip(dot(a, o, y, col, size=7.5), f"{g} — {lab}: TEST OR {o:.2f} (jackknife 95% CI {lo:.2f}–{hi:.2f})")
            rel = "" if o == base else f"  ({(o / base - 1) * 100:+.0f}%)"
            a.text(hi * 1.12, y, f"{fmt_or(o)}{rel}", va="center", fontsize=8.3,
                   color=INK if o == base else INK2, fontweight="bold" if o == base else "normal")
            yt.append(y)
            ylab.append(lab)
            y += 1
        y += 0.9
    a.set_yticks(yt, ylab, fontsize=8.2)
    a.tick_params(axis="y", length=0, labelcolor=INK)
    a.set_ylim(y - 0.6, -1.35)
    log_or_axis(a, [0.5, 1, 2, 5, 10, 20, 50])
    a.set_xlim(0.4, 60)
    a.axvline(1, color=INK2, lw=1, zorder=1)
    style_axis(a, "x")
    a.spines["left"].set_visible(False)
    a.set_xlabel("TEST odds ratio, top 0.1% (log scale; bars = jackknife 95% CI)")
    a.set_title("a  Shuffle the targets: OR barely moves. Permute CN: OR falls (CN-only: n.s.)", fontsize=9.5)

    # (b) per-regulator AUROC
    rows = [("g1_t100_w0__a3_r095_s2_w0", "CN-only: real", BLUE),
            ("g1_t100_w0__a3_r095_s2_w0_C1", "CN-only: C1 shuffle", AQUA),
            ("g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2", "context layer: real", BLUE),
            ("g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1", "context layer: C1 shuffle", AQUA)]
    yy = [0, 1, 2.6, 3.6]
    for (run, lab, col), y0 in zip(rows, yy):
        au = final_json(run)["auroc"]
        b.plot([0.5, au["mean_auroc"]], [y0, y0], color=GRID, lw=2, zorder=1, solid_capstyle="butt")
        tip(dot(b, au["mean_auroc"], y0, col, size=7.5),
            f"{lab}: mean per-regulator AUROC {au['mean_auroc']:.3f} over {au['n_regulators']} regulators")
        b.text(au["mean_auroc"] + 0.012, y0, f"{au['mean_auroc']:.3f}", va="center", fontsize=8.3, color=INK)
    b.set_yticks(yy, [r[1] for r in rows], fontsize=8.2)
    b.tick_params(axis="y", length=0, labelcolor=INK)
    b.set_ylim(4.3, -0.9)
    b.set_xlim(0.45, 0.72)
    b.axvline(0.5, color=INK2, lw=1, zorder=1)
    b.text(0.503, -0.6, "chance", fontsize=7.8, color=INK2)
    b.set_xticks([0.45, 0.5, 0.55, 0.6, 0.65, 0.7])
    style_axis(b, "x")
    b.spines["left"].set_visible(False)
    n_reg = final_json("g1_t100_w0__a3_r095_s2_w0")["auroc"]["n_regulators"]
    b.set_xlabel(f"mean per-regulator AUROC ({n_reg} regulators)")
    b.set_title("b  Within a regulator, CN ranks targets at chance", fontsize=9.5)
    b.text(0.507, 3.1, "the context layer's 0.65 comes from target\nexpression and survives the shuffle",
           fontsize=7.6, color=INK2, va="center", bbox=dict(fc=SURFACE, ec="none", pad=1.5))
    cmh = final_json("g1_t100_w0__a3_r095_s2_w0")["cmh_regulator_x_exprdecile"]
    fig.text(0.015, 0.93, "What signal there is: which regulators, not which targets (HCT116 TEST)", fontsize=10.5,
             fontweight="bold", color=INK)
    # paired Δlog-OR values from results/final/summary.md ("paired" table)
    fig.text(0.015, 0.018, f"C1 shuffles each regulator's scores among its targets (within expression bins); NC1/C2 "
             f"permute CN. CMH OR (CN-only) = {cmh['cmh_or']:.2f} (95% CI {cmh['cmh_ci95'][0]:.2f}–"
             f"{cmh['cmh_ci95'][1]:.2f}).\nReal vs permuted CN, paired Δlog-OR: CN-only +1.10 (CI −0.66 to 2.85, "
             f"not significant); context layer +0.55 (CI 0.07–1.03).",
             fontsize=7.5, color=MUTED)
    save(fig, "what_the_signal_is", "What signal there is concerns which regulators, not which targets",
         "Shuffling scores among each regulator's targets leaves the TEST odds ratio almost unchanged (3.85 → 3.76; "
         "context layer 18.4 → 17.3) and per-regulator AUROC at 0.50, whereas permuting copy number, which changes "
         "which regulators score high, lowers it (1.28; context fusion 10.6). For the CN-only tool that drop is not "
         "significant (paired Δlog-OR +1.10, CI −0.66 to 2.85); for the context layer it is (+0.55, CI 0.07–1.03).",
         ["results/final/g1_t100_w0__a3_r095_s2_w0/hct116_TEST.json",
          "results/final/g1_t100_w0__a3_r095_s2_w0_C1/hct116_TEST.json",
          "results/final/g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2/hct116_TEST.json",
          "results/final/g1_t100_w0__a3_r095_s2_w0__i1_fuse_g2_C1/hct116_TEST.json",
          "results/final/NC1/hct116_TEST.json", "results/final/C2_g2/hct116_TEST.json"])


# ================================================================= F8 journey
MILESTONES = [  # (commit, label, level: + above / - below the line; |2| = second tier)
    ("e31c6f7", "setup", 1),
    ("aacc223", "literature + data\nprofile", -1),
    ("a6ccaeb", "design frozen after\n3 review rounds", 1),
    ("1a4a211", "harness + model\nbuilt", -1),
    ("57f2de3", "TEST labels\nsealed", 2),
    ("f2deb85", "M1 weak (OR 1.6)\nvs hub baselines", -2),
    ("6288007", "G1 cis gate\nadopted (3.5)", 1),
    ("c968b1d", "A3 attribution\nadopted (5.0)", -1),
    ("d439c58", "context fusion fails\npair-shuffle control", 1),
    ("e9085c7", "single TEST run\n(LOCKed)", -1),
    ("b0827fd", "report verified\n(4 errors fixed)", 1),
]
PHASES = [("e31c6f7", "a6ccaeb", "plan + design"), ("74e9e6c", "57f2de3", "build"),
          ("f2deb85", "d439c58", "tune on TRAIN"), ("2b4242c", "c27a904", "final run + report")]


def fig_journey():
    out = subprocess.run(["git", "-C", str(ROOT), "log", "--date=iso", "--pretty=%ad|%h|%s"],
                         capture_output=True, text=True, check=True).stdout.strip().splitlines()
    commits = {}
    for ln in out:
        d, h, s = ln.split("|", 2)
        commits[h] = (datetime.strptime(d.strip(), "%Y-%m-%d %H:%M:%S %z"), s)
    t0 = commits["e31c6f7"][0]
    end = commits["c27a904"][0]
    hrs = lambda h: (commits[h][0] - t0).total_seconds() / 3600  # noqa: E731
    fig, ax = plt.subplots(figsize=(11, 4.0))
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.15)
    span = (end - t0).total_seconds() / 3600
    ax.plot([0, span], [0, 0], color=AXIS, lw=1.2, zorder=1)
    # all commits as ticks
    for h, (d, s) in commits.items():
        x = (d - t0).total_seconds() / 3600
        if 0 <= x <= span:
            tk, = ax.plot([x, x], [-0.07, 0.07], color=MUTED, lw=0.8, zorder=2)
            tip(tk, f"{d:%H:%M} {h}: {s}")
    # phases as a thin band below
    for a_, b_, lab in PHASES:
        x0, x1 = hrs(a_), hrs(b_)
        ax.add_patch(Rectangle((x0, -1.30), x1 - x0, 0.08, fc=BLUE_150, ec="none", zorder=1))
        ax.text((x0 + x1) / 2, -1.37, lab, ha="center", va="top", fontsize=7.8, color=INK2)
    for h, lab, lvl in MILESTONES:
        x = hrs(h)
        up = lvl > 0
        yl = {1: 0.38, 2: 0.78, -1: -0.30, -2: -0.70}[lvl]
        ax.plot([x, x], [0, yl - (0.04 if up else -0.04)], color=AXIS, lw=0.8, zorder=1)
        col = ORANGE if h in ("57f2de3", "e9085c7") else BLUE
        tip(dot(ax, x, 0, col, size=7.5, z=5), f"{commits[h][0]:%Y-%m-%d %H:%M} UTC · {h}: {commits[h][1]}")
        ax.text(x, yl, lab, ha="center", va="bottom" if up else "top", fontsize=7.8, color=INK, linespacing=1.2,
                bbox=dict(fc=SURFACE, ec="none", pad=0.5))
    ticks = list(range(0, int(span) + 1))
    ax.set_xticks(ticks, [f"{(t0.hour + t) % 24:02d}:00" if False else f"+{t} h" for t in ticks])
    ax.set_xlim(-0.25, span + 0.25)
    ax.set_ylim(-1.65, 1.25)
    ax.set_yticks([])
    for s_ in ("top", "right", "left"):
        ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.set_xlabel(f"hours since project start ({t0:%Y-%m-%d %H:%M} UTC); ticks on the line = commits")
    fig.text(0.03, 0.93, "Project journey, placed by git commit time", fontsize=10.5, fontweight="bold", color=INK)
    lx = 0.62
    for i, (c, t) in enumerate([(BLUE, "milestone"), (ORANGE, "TEST seal / single TEST run")]):
        fig.lines.append(plt.Line2D([lx + i * 0.11], [0.935], marker="o", ms=6.5, mfc=c, mec=SURFACE,
                                    transform=fig.transFigure))
        fig.text(lx + i * 0.11 + 0.01, 0.935, t, fontsize=8, color=INK2, va="center")
    save(fig, "journey", "Project journey",
         f"The whole project ran in about {span:.1f} hours of commits on 27 Oct 2026: design frozen after three "
         "review rounds, TEST labels sealed before any evaluation, two tuning steps adopted on TRAIN, a context "
         "layer rejected by its control, then one locked TEST run and an independent report check.",
         ["git log --date=iso"])


def main():
    fig_pipeline()
    fig_benchmark_structure()
    fig_split()
    fig_tuning_path()
    fig_test_results()
    fig_sign_concordance()
    fig_what_the_signal_is()
    fig_journey()
    (OUT / "figures.json").write_text(json.dumps(
        {"generated_by": "scripts/figures/make_figures.py", "theme": "light",
         "palette": {"categorical": CATEGORICAL, "ink": INK, "secondary_ink": INK2, "muted": MUTED,
                     "grid": GRID, "surface": SURFACE},
         "figures": MANIFEST}, indent=1, ensure_ascii=False))
    print("wrote figures.json")


if __name__ == "__main__":
    main()
