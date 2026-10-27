"""Regulator-level, sign, incremental within-regulator value of |t| over target priors; focal (diagnostic only)."""
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
T, scis, cis = load('g1_t100_w0__a3_r095_s2_w0', 't')
S, _, _ = load('g1_t100_w0__a3_r095_s2_w0', 'score')
A, _, _ = load('g1_t100_w0__a3_r095_s2_w0', 'alpha')
F, _, _ = load('focal', 't'); B6, _, _ = load('B6'); B3, _, _ = load('B3'); B4, _, _ = load('B4')
ev = bm.Evaluation(regs, tgts, pos, S, scis); nt, nr = ev.nt, ev.nr
trans = ev.key > -2
P = np.zeros((nr, nt), bool); pl = ev.pos_lin(); P.ravel()[pl] = True
eh = tgts['expr_hct116'].values; deg = P.sum(1)
g1 = (cis > 0) & (cis**2 >= 100)
ctx = regs['expr_context'].values; own = regs['own_class'].values
sel = ev.topk(0.001); nsel = np.bincount(sel // nt, minlength=nr)
q999 = np.array([np.nanpercentile(np.abs(S[i, trans[i]]), 99.9) for i in range(nr)])
print('n TRAIN regs', nr, 'G1 pass', g1.sum(), 'G1 & ctx>=1', (g1 & (ctx >= 1)).sum(), 'G1 & ctx<1', (g1 & (ctx < 1)).sum())
print('top-K pairs from G1 regs with ctx<1: %.3f ; with ctx<2: %.3f' % (nsel[g1 & (ctx < 1)].sum()/len(sel), nsel[g1 & (ctx < 2)].sum()/len(sel)))
print('hit rate (positives per universe pair) by ctx expression of regulator:')
for lo, hi in [(-1, 1), (1, 3), (3, 5), (5, 7), (7, 20)]:
    m = (ctx >= lo) & (ctx < hi); print(f'  ctx [{lo},{hi}) regs {m.sum():5d} rate {deg[m].sum()/ev.n_i[m].sum():.5f}  (G1 subset rate {deg[m&g1].sum()/max(1,ev.n_i[m&g1].sum()):.5f}, topK share {nsel[m].sum()/len(sel):.3f})')
print('degree by own class:', pd.Series(deg).groupby(own).mean().round(1).to_dict())
print('Spearman(regulator degree, x) all TRAIN regs / G1-pass regs:')
b3 = np.abs(B3[:, 0])
for name, x in [('cis_t', cis), ('B3=|cis t|*SD(x)', b3), ('expr_context', ctx), ('q99.9|S| row', q999), ('n topK', nsel), ('n_lines|log2CN|>=1', regs['cn_stat'].values)]:
    print(f'  {name:22s} {spearmanr(x, deg)[0]:.3f} / {spearmanr(x[g1], deg[g1])[0]:.3f}')
# ---------- sign ----------
st = sel % nt; sr = sel // nt
tv = T.ravel()[sel]; hit = P.ravel()[sel]
for nm, m in [('t>0 (pred KD down)', tv > 0), ('t<0 (pred KD up)', tv < 0)]:
    print(f'top-K {nm}: n {m.sum()} hit rate {hit[m].mean():.4f}')
# sign of positives overall among G1 regs: does sign(t) agree w/ LFC beyond top-K? within top 1% of each G1 reg
# ---------- incremental within-reg AUROC ----------
def auc(x, y):
    r = rankdata(x); n1 = y.sum(); n0 = len(y) - n1
    return (r[y].sum() - n1*(n1+1)/2) / (n1*n0)
rows = np.flatnonzero((deg >= 10) & g1)
At = np.where(trans, np.abs(T), np.nan); tgtresp = np.nanmean(At, 0)
sdB = B4[0]
out = {}
for i in rows:
    m = trans[i] & np.isfinite(T[i]); y = P[i, m]
    rt = rankdata(np.abs(T[i, m])) / m.sum(); re = rankdata(eh[m]) / m.sum()
    cand = {'|t|': rt, 'eh': re, 'eh+0.25|t|': re + 0.25*rt, 'eh+0.5|t|': re + 0.5*rt, 'eh*|t| (rank prod)': re*rt,
            '|t|*SD_B': np.abs(T[i, m])*sdB[m], '|r|*sqrt(TPM)': np.abs(T[i, m])*np.sqrt(2**eh[m]-1),
            '|S|(alpha)': np.abs(S[i, m]), 'alpha': A[i, m],
            '|t_focal|': np.abs(np.nan_to_num(F[i, m])), 'eh+0.5|t_focal|': re + 0.5*rankdata(np.abs(np.nan_to_num(F[i, m])))/m.sum(),
            'tgtresp': tgtresp[m], 'eh+0.5tgtresp': re + 0.5*rankdata(tgtresp[m])/m.sum(),
            't signed (+ high)': T[i, m]}
    for k, v in cand.items(): out.setdefault(k, []).append(auc(v, y))
print(f'within-regulator AUROC, G1-pass regs with >=10 positives (n={len(rows)}):')
for k, v in out.items():
    v = np.array(v); print(f'  {k:22s} mean {v.mean():.4f}  median {np.median(v):.4f}  frac>eh {np.mean(v > np.array(out["eh"])):.2f}')
# focal vs LOCO: corr with B6, within G1
cf = [np.corrcoef(np.nan_to_num(F[i, trans[i]]), np.nan_to_num(B6[i, trans[i]]))[0, 1] for i in np.flatnonzero(g1)[:300]]
cl = [np.corrcoef(np.nan_to_num(F[i, trans[i]]), np.nan_to_num(T[i, trans[i]]))[0, 1] for i in np.flatnonzero(g1)[:300]]
print('focal t row: median corr with B6 %.3f, with LOCO t %.3f' % (np.nanmedian(cf), np.nanmedian(cl)))
print('alpha: frac of G1 trans pairs with alpha<0.5: %.3f; among topK %.3f' % ((A[g1][trans[g1]] < 0.5).mean(), (A.ravel()[sel] < 0.5).mean()))
# segment sizes proxy: number of TRAIN regulators sharing near-identical t rows? skip
# rows with large |t| magnitudes: max |t| by arm
print('median of row-max |t| (G1): %.1f; t* call threshold 3.0; fraction of G1 trans pairs with |t|>3: %.3f' % (np.nanmedian(np.nanmax(At[g1], 1)), (At[g1] > 3).mean()))
