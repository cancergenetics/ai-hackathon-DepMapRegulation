"""Within-regulator AUROC of |t| by regulator degree / instrument strength (diagnostic; HCT116 TRAIN)."""
import numpy as np, pandas as pd, sys
from scipy.stats import rankdata, spearmanr
sys.path.insert(0, '/Users/colmclaude/ClaudeExperiments/DepMapTR/src')
from depmaptr import benchmark as bm
R = '/Users/colmclaude/ClaudeExperiments/DepMapTR/data/results/scores/'
bench = bm.load_bench('hct116', 'TRAIN'); regs, tgts, pos = bench.regs, bench.tgts, bench.pos
z = np.load(R + 'g1_t100_w0__a3_r095_s2_w0/scores.npz')
sc = {'regulators': z['regulators'].astype(str), 'targets': z['targets'].astype(str), 's_cis': z['s_cis'], 't': z['t']}
T, scis, _, _ = bm.align_scores(sc, regs['regulator'], tgts['target'], 't')
cis = z['cis_t'][pd.Index(sc['regulators']).get_indexer(regs['regulator'])]
ev = bm.Evaluation(regs, tgts, pos, T, scis); nt, nr = ev.nt, ev.nr; trans = ev.key > -2
P = np.zeros((nr, nt), bool); P.ravel()[ev.pos_lin()] = True; deg = P.sum(1); eh = tgts['expr_hct116'].values
def auc(x, y):
    r = rankdata(x); n1 = y.sum(); n0 = len(y) - n1
    return (r[y].sum() - n1*(n1+1)/2) / (n1*n0)
rows = np.flatnonzero(deg >= 10)
A = np.array([auc(np.abs(np.nan_to_num(T[i, trans[i]])), P[i, trans[i]]) for i in rows])
# eh-adjusted: AUROC within HCT116 expression deciles (rank |t| within decile)
d = pd.qcut(rankdata(eh), 10, labels=False)
def adj(i):
    m = trans[i]; x = np.abs(np.nan_to_num(T[i, m])); dd = d[m]; rr = np.empty(m.sum())
    for k in range(10):
        w = dd == k; rr[w] = rankdata(x[w]) / w.sum()
    return auc(rr, P[i, m])
Aa = np.array([adj(i) for i in rows])
g1 = (cis > 0) & (cis**2 >= 100)
for lo, hi in [(10, 30), (30, 100), (100, 300), (300, 10**6)]:
    m = (deg[rows] >= lo) & (deg[rows] < hi)
    mg = m & g1[rows]
    print(f'degree [{lo},{hi}): n {m.sum():3d} AUROC|t| {A[m].mean():.4f} expr-adj {Aa[m].mean():.4f} | G1 n {mg.sum():3d} AUROC {A[mg].mean() if mg.any() else np.nan:.4f} adj {Aa[mg].mean() if mg.any() else np.nan:.4f}')
print('spearman(AUROC, cis t) over regs deg>=10: %.3f; G1 only %.3f' % (spearmanr(A, cis[rows])[0], spearmanr(A[g1[rows]], cis[rows][g1[rows]])[0]))
se = np.sqrt(np.mean([(a*(1-a)) for a in A]) / 30)
print('mean AUROC all %.4f (se approx via regulators sd/sqrt(n)): %.4f; frac>0.55 %.2f; frac<0.45 %.2f' % (A.mean(), A.std()/np.sqrt(len(A)), (A > 0.55).mean(), (A < 0.45).mean()))
