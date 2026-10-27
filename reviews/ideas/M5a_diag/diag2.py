"""Within-regulator diagnostics, arm signatures, alternative CMH strata (diagnostic only; HCT116 TRAIN)."""
import numpy as np, pandas as pd, sys
from scipy.stats import spearmanr, rankdata
sys.path.insert(0, '/Users/colmclaude/ClaudeExperiments/DepMapTR/src')
from depmaptr import benchmark as bm
R = '/Users/colmclaude/ClaudeExperiments/DepMapTR/data/results/scores/'
bench = bm.load_bench('hct116', 'TRAIN'); regs, tgts, pos = bench.regs, bench.tgts, bench.pos
def load(run, key='score'):
    z = np.load(R + run + '/scores.npz')
    sc = {'regulators': z['regulators'].astype(str), 'targets': z['targets'].astype(str), key: z[key], 's_cis': z['s_cis']}
    S, s, _, _ = bm.align_scores(sc, regs['regulator'], tgts['target'], key)
    ci = pd.Index(sc['regulators']).get_indexer(regs['regulator'])
    return S, s, (z['cis_t'][ci] if 'cis_t' in z.files else None)
T, scis, cis = load('g1_t100_w0__a3_r095_s2_w0', 't')     # raw reduced-form t
S, _, _ = load('g1_t100_w0__a3_r095_s2_w0', 'score')
B2, _, _ = load('B2'); B6, _, _ = load('B6')
ev = bm.Evaluation(regs, tgts, pos, S, scis)
nt, nr = ev.nt, ev.nr
trans = ev.key > -2
P = np.zeros((nr, nt), bool); pl = ev.pos_lin(); P.ravel()[pl] = True
eh = tgts['expr_hct116'].values
npos_r = P.sum(1)
# ---------- within-regulator AUROC (regs with >=10 positives) ----------
def auc(x, y):
    r = rankdata(x); n1 = y.sum(); n0 = len(y) - n1
    return (r[y].sum() - n1*(n1+1)/2) / (n1*n0)
rows = np.flatnonzero(npos_r >= 10)
g1 = (cis > 0) & (cis**2 >= 100)
res = {k: [] for k in ['absT', 'absT_g1', 'eh', 'absB2', 'absB6', 'absT_resid_eh', 'tgtresp']}
# target CN-responsiveness: mean |t| over all TRAIN regulators (trans entries)
At = np.where(trans, np.abs(T), np.nan)
tgtresp = np.nanmean(At, 0)
for i in rows:
    m = trans[i] & np.isfinite(T[i]); y = P[i, m]
    res['absT'].append(auc(np.abs(T[i, m]), y))
    if g1[i]: res['absT_g1'].append(auc(np.abs(T[i, m]), y))
    res['eh'].append(auc(eh[m], y))
    b = B2[i, m]; 
    if np.isfinite(b).all(): res['absB2'].append(auc(np.abs(b), y))
    b6 = B6[i, m]
    if np.isfinite(b6).all(): res['absB6'].append(auc(np.abs(b6), y))
    # |t| within HCT116-expression decile strata: AUROC of |t| rank-normalised within decile
    d = pd.qcut(eh[m], 10, labels=False, duplicates='drop')
    rr = np.zeros(m.sum())
    for k in np.unique(d):
        w = d == k; rr[w] = rankdata(np.abs(T[i, m][w])) / w.sum()
    res['absT_resid_eh'].append(auc(rr, y))
    res['tgtresp'].append(auc(tgtresp[m], y))
print('within-regulator AUROC (mean over regs with >=10 pos; n=%d; G1-pass n=%d)' % (len(rows), len(res['absT_g1'])))
for k, v in res.items(): print(f'  {k:15s} {np.mean(v):.4f}  (n={len(v)}, frac>0.5={np.mean(np.array(v)>0.5):.2f})')
# ---------- correlation of within-reg |t| with target expression & responsiveness ----------
rs_eh, rs_resp = [], []
rng = np.random.default_rng(0)
for i in rng.choice(np.flatnonzero(g1), 200, replace=False):
    m = trans[i] & np.isfinite(T[i])
    rs_eh.append(spearmanr(np.abs(T[i, m]), eh[m])[0]); rs_resp.append(spearmanr(np.abs(T[i, m]), tgtresp[m])[0])
print('median spearman(|t| row, target HCT116 expr) over 200 G1 regs:', np.median(rs_eh).round(3))
print('median spearman(|t| row, target mean|t| across regs):', np.median(rs_resp).round(3))
print('spearman(target mean|t|, target benchmark degree):', spearmanr(tgtresp, P.sum(0))[0].round(3))
print('spearman(target mean|t|, expr_hct116):', spearmanr(tgtresp, eh)[0].round(3), ' vs SD across S2 (B4):')
B4, _, _ = load('B4'); sd_t = B4[0]
print('spearman(target mean|t|, SD_S2):', spearmanr(tgtresp, sd_t)[0].round(3))
# ---------- arm signatures: corr of t row with B6 arm row, and pairwise among same-arm G1 regs ----------
arm = regs['arm'].values
cr = []
for i in np.flatnonzero(g1):
    m = trans[i] & np.isfinite(T[i]) & np.isfinite(B6[i])
    if m.sum() > 1000: cr.append(np.corrcoef(T[i, m], B6[i, m])[0, 1])
cr = np.array(cr); print('G1 regs: corr(t row, B6 arm row): median %.3f, q10 %.3f, q90 %.3f, frac>0.7 %.2f' % (np.median(cr), *np.percentile(cr, [10, 90]), (cr > 0.7).mean()))
# top-K composition by arm & distinct segments
sel = ev.topk(0.001); sr = sel // nt
for a in ['1q', '1p', '10q', '12p', '13q']:
    ri = np.flatnonzero((arm == a))
    rsel = np.intersect1d(ri, np.unique(sr))
    # overlap of top-K target sets among regs of this arm
    sets = [set((sel[sr == i] % nt).tolist()) for i in rsel]
    allt = set().union(*sets) if sets else set()
    tot = sum(len(s) for s in sets)
    print(f'arm {a}: {len(rsel)} regs in topK, {tot} pairs, {len(allt)} distinct targets (redundancy {tot/max(1,len(allt)):.2f}); arm positives rate {P[ri][trans[ri]].mean():.5f}')
# ---------- alternative CMH: regulator x HCT116-expression decile (diagnostic of pair-specific signal) ----------
def cmh_alt(Smat, label):
    e2 = bm.Evaluation(regs, tgts.assign(expr_decile=pd.qcut(rankdata(eh), 10, labels=False) + 1), pos, Smat, scis)
    s = e2.topk(0.001); pl2 = e2.pos_lin()
    c1 = ev.cmh(ev.topk(0.001) if Smat is S else bm.Evaluation(regs, tgts, pos, Smat, scis).topk(0.001), pl, 'regulator')['cmh_or'] if False else None
    print(f'  {label:10s} CMH(reg x HCT116-expr decile) = {e2.cmh(s, pl2)["cmh_or"]:.3f}')
print('Alternative CMH strata (diagnostic only):')
cmh_alt(S, 'incumbent'); cmh_alt(B2, 'B2')
B5, _, _ = load('B5'); cmh_alt(B5, 'B5')
B4e, _, _ = load('B4e'); cmh_alt(B4e, 'B4e')
