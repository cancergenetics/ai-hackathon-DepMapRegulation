"""TRAIN-only diagnostics for M5a ideas review (no tuning; HCT116 TRAIN)."""
import numpy as np, pandas as pd, sys
from scipy.stats import spearmanr, rankdata
sys.path.insert(0, '/Users/colmclaude/ClaudeExperiments/DepMapTR/src')
from depmaptr import benchmark as bm

R = '/Users/colmclaude/ClaudeExperiments/DepMapTR/data/results/scores/'
bench = bm.load_bench('hct116', 'TRAIN')
regs, tgts, pos = bench.regs, bench.tgts, bench.pos

def load(run, key='score'):
    z = np.load(R + run + '/scores.npz')
    sc = {'regulators': z['regulators'].astype(str), 'targets': z['targets'].astype(str), key: z[key], 's_cis': z['s_cis']}
    S, s, _, _ = bm.align_scores(sc, regs['regulator'], tgts['target'], key)
    extra = {}
    if 'cis_t' in z.files:
        ci = pd.Index(sc['regulators']).get_indexer(regs['regulator']); extra['cis_t'] = z['cis_t'][ci]
    return S, s, extra

S, scis, ex = load('g1_t100_w0__a3_r095_s2_w0')
ev = bm.Evaluation(regs, tgts, pos, S, scis)
sel = ev.topk(0.001); pl = ev.pos_lin()
nt = ev.nt
sr, st = sel // nt, sel % nt
print('K', len(sel), 'distinct regs', len(np.unique(sr)), 'distinct targets', len(np.unique(st)))
rc = np.bincount(sr, minlength=ev.nr)
o = np.argsort(-rc)
print('top10 reg share', rc[o[:10]].sum()/len(sel), 'top30', rc[o[:30]].sum()/len(sel), 'top100', rc[o[:100]].sum()/len(sel))
arm = regs['arm'].values
ac = pd.Series(arm[sr]).value_counts()
print('top-K by regulator arm (share):'); print((ac/len(sel)).head(8).round(3).to_string())
print('n arms in topK', len(ac))
# positives per top regulators
s_, p_, a_, both = ev.per_reg(sel, pl)
top = o[:15]
df = pd.DataFrame({'reg': regs['regulator'].values[top], 'arm': arm[top], 'nsel': s_[top], 'npos': p_[top], 'hit': a_[top],
                   'own': regs['own_class'].values[top], 'cis_t': ex['cis_t'][top], 'expr_ctx': regs['expr_context'].values[top]})
df['prec'] = df.hit/df.nsel; df['base'] = df.npos/ev.n_i[top]
print(df.round(3).to_string())
# overall precision vs base
print('pooled precision in topK', len(both)/len(sel), 'base rate', len(pl)/ev.N)
# precision expected if within-regulator selection were random: sum_r nsel_r * base_r / K
base_r = p_/ev.n_i
print('expected precision under within-reg random', (s_*base_r).sum()/len(sel))
# and if additionally target expression mattered: strat by hct116 expr decile within reg
eh = tgts['expr_hct116'].values
print('\nTarget HCT116 expression: topK targets vs universe vs positives')
q = lambda v: np.percentile(v, [10,25,50,75,90]).round(2)
print(' universe', q(eh), 'frac<1:', (eh<1).mean().round(3))
print(' topK   ', q(eh[st]), 'frac<1:', (eh[st]<1).mean().round(3))
print(' pos    ', q(eh[pl % nt]), 'frac<1:', (eh[pl % nt]<1).mean().round(3))
# positive rate by HCT116 expression bin (pairs from all TRAIN regs)
bins = np.array([-1, 0.5, 1, 2, 3, 4, 5, 6, 7, 20])
tb = np.digitize(eh, bins) - 1
npos_t = np.bincount(pl % nt, minlength=nt)
ntb = np.bincount(tb, minlength=len(bins)); posb = np.bincount(tb, weights=npos_t, minlength=len(bins))
selb = np.bincount(tb[st], minlength=len(bins))
print(' bin  ntargets  pos_rate_per_target  share_of_topK')
for i in range(len(bins)-1):
    if ntb[i]: print(f'  [{bins[i]},{bins[i+1]}) {ntb[i]:6d} {posb[i]/ntb[i]:8.2f} {selb[i]/len(sel):.3f}')
# target concentration
tcnt = np.bincount(st, minlength=nt)
ot = np.argsort(-tcnt)
print('top-50 targets share of topK', tcnt[ot[:50]].sum()/len(sel), 'top-200', tcnt[ot[:200]].sum()/len(sel))
print('spearman(target topK count, target benchmark degree) over targets with >=1 topK:',
      spearmanr(tcnt[tcnt>0], npos_t[tcnt>0])[0])
print('spearman(target benchmark degree, expr_hct116) all targets', spearmanr(npos_t, eh)[0])
print('spearman(target benchmark degree, expr_mean_s2) all targets', spearmanr(npos_t, tgts['expr_mean_s2'].values)[0])
np.save('/private/tmp/claude-503/-Users-colmclaude-ClaudeExperiments-DepMapTR/130bbdd0-af4a-481e-a7c7-41816a9e851c/scratchpad/npos_t.npy', npos_t)
