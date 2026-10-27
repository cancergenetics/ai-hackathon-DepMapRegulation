"""Sign asymmetry within regulator, segment sizes (diagnostic only; HCT116 TRAIN)."""
import numpy as np, pandas as pd, sys
sys.path.insert(0, '/Users/colmclaude/ClaudeExperiments/DepMapTR/src')
from depmaptr import benchmark as bm
R = '/Users/colmclaude/ClaudeExperiments/DepMapTR/data/results/scores/'
bench = bm.load_bench('hct116', 'TRAIN'); regs, tgts, pos = bench.regs, bench.tgts, bench.pos
z = np.load(R + 'g1_t100_w0__a3_r095_s2_w0/scores.npz')
sc = {k: z[k] for k in ['regulators', 'targets', 's_cis']}; sc['regulators'] = sc['regulators'].astype(str); sc['targets'] = sc['targets'].astype(str)
sc['score'] = z['score']; S, scis, _, _ = bm.align_scores(sc, regs['regulator'], tgts['target'])
sc['t'] = z['t']; T, _, _, _ = bm.align_scores(sc, regs['regulator'], tgts['target'], 't')
ev = bm.Evaluation(regs, tgts, pos, S, scis); nt, nr = ev.nt, ev.nr
P = np.zeros((nr, nt), bool); P.ravel()[ev.pos_lin()] = True
sel = ev.topk(0.001); sr, st = sel // nt, sel % nt
tv = T.ravel()[sel]; hit = P.ravel()[sel]; eh = tgts['expr_hct116'].values
# MH OR of hit for t>0 vs t<0 stratified by regulator (top-K pairs only)
num = den = 0.0
for r in np.unique(sr):
    m = sr == r; a = (hit & (tv > 0))[m].sum(); b = ((~hit) & (tv > 0))[m].sum(); c = (hit & (tv < 0))[m].sum(); d = ((~hit) & (tv < 0))[m].sum()
    n = a + b + c + d
    if n: num += a * d / n; den += b * c / n
print('within-regulator MH OR (hit | t>0 vs t<0) among top-K: %.3f' % (num / den))
print('median HCT116 expr of top-K targets: t>0 %.2f, t<0 %.2f' % (np.median(eh[st[tv > 0]]), np.median(eh[st[tv < 0]])))
print('fraction of top-K regulators whose top-K pairs are >90%% one sign: %.2f' % np.mean([max((tv[sr == r] > 0).mean(), (tv[sr == r] < 0).mean()) > 0.9 for r in np.unique(sr)]))
# all G1 regs, all trans pairs: hit rate by sign of t in |t|>=3
cis = z['cis_t'][pd.Index(sc['regulators']).get_indexer(regs['regulator'])]
g1 = (cis > 0) & (cis**2 >= 100)
trans = ev.key > -2
for nm, m in [('t>=3', (T >= 3)), ('t<=-3', (T <= -3)), ('|t|<3', (np.abs(T) < 3))]:
    mm = m & trans & g1[:, None]
    print(f'G1 regs, pairs {nm}: n {mm.sum()} hit rate {P[mm].mean():.5f}')
