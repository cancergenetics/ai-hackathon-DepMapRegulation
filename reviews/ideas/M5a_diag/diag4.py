"""Focal instrument, arm allocation, eh-matched tail precision (diagnostic only; HCT116 TRAIN)."""
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
T, scis, cis = load('g1_t100_w0__a3_r095_s2_w0', 't'); S, _, _ = load('g1_t100_w0__a3_r095_s2_w0')
F, _, cisF = load('focal', 't')
ev = bm.Evaluation(regs, tgts, pos, S, scis); nt, nr = ev.nt, ev.nr; trans = ev.key > -2
P = np.zeros((nr, nt), bool); pl = ev.pos_lin(); P.ravel()[pl] = True
eh = tgts['expr_hct116'].values; deg = P.sum(1); g1 = (cis > 0) & (cis**2 >= 100); ctx = regs['expr_context'].values
print('focal cis t: median %.2f (LOCO %.2f); n F>=100 focal %d, LOCO %d; both %d' % (np.median(cisF), np.median(cis), (cisF**2*(cisF>0) >= 100).sum(), g1.sum(), ((cisF**2*(cisF>0) >= 100) & g1).sum()))
print('Spearman(degree, focal cis t) all %.3f; among G1 %.3f' % (spearmanr(cisF, deg)[0], spearmanr(cisF[g1], deg[g1])[0]))
def auc(x, y):
    r = rankdata(x); n1 = y.sum(); n0 = len(y) - n1
    return (r[y].sum() - n1*(n1+1)/2) / (n1*n0)
for lo, hi in [(0, 5), (5, 10), (10, 100)]:
    m = (cisF >= lo) & (cisF < hi) & (deg >= 10) & g1
    a_f = [auc(np.abs(np.nan_to_num(F[i, trans[i]])), P[i, trans[i]]) for i in np.flatnonzero(m)]
    a_l = [auc(np.abs(np.nan_to_num(T[i, trans[i]])), P[i, trans[i]]) for i in np.flatnonzero(m)]
    print(f'G1 regs with focal cis t in [{lo},{hi}) & deg>=10: n={m.sum()}, within-reg AUROC focal |t| {np.mean(a_f) if a_f else np.nan:.4f}, LOCO |t| {np.mean(a_l) if a_l else np.nan:.4f}')
# arm allocation vs arm hit rate
sel = ev.topk(0.001); sr = sel // nt
arm = regs['arm'].values; ua = np.unique(arm)
share = np.array([np.isin(sr, np.flatnonzero(arm == a)).sum() / len(sel) for a in ua])
rate = np.array([deg[arm == a].sum() / ev.n_i[arm == a].sum() for a in ua])
rateg1 = np.array([deg[(arm == a) & g1].sum() / max(1, ev.n_i[(arm == a) & g1].sum()) for a in ua])
print('arms: Spearman(topK share, arm hit rate) %.3f; (topK share, G1 arm hit rate) %.3f' % (spearmanr(share, rate)[0], spearmanr(share, rateg1)[0]))
o = np.argsort(-share)
print(pd.DataFrame({'arm': ua[o], 'share': share[o].round(3), 'rate_all': rate[o].round(5), 'rate_G1': rateg1[o].round(5)}).head(10).to_string(index=False))
# eh-matched expected precision for top-K: per pair, base rate of its (regulator x HCT116-expr-decile) stratum
d = pd.qcut(rankdata(eh), 10, labels=False)
st = sel % nt
num = np.zeros((nr, 10)); den = np.zeros((nr, 10))
for k in range(10):
    cols = d == k
    num[:, k] = (P[:, cols] & trans[:, cols]).sum(1); den[:, k] = trans[:, cols].sum(1)
br = num / np.maximum(den, 1)
exp_prec = br[sr, d[st]].mean(); obs = P.ravel()[sel].mean()
exp_reg = (deg / ev.n_i)[sr].mean()
print(f'top-K precision: observed {obs:.5f}; expected if within-reg random {exp_reg:.5f}; expected if within-reg random but HCT116-expr-decile matched {exp_prec:.5f}')
# what B4e-within-reg would give with incumbent's regulator allocation: precision of top-n_r targets by eh for each reg
nsel = np.bincount(sr, minlength=nr)
hits = 0
for i in np.flatnonzero(nsel):
    m = np.flatnonzero(trans[i]); top = m[np.argsort(-eh[m])[:nsel[i]]]; hits += P[i, top].sum()
print(f'precision with incumbent regulator allocation but targets ranked by HCT116 expression: {hits/len(sel):.5f}')
