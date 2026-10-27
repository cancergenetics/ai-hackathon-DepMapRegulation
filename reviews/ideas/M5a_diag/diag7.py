"""Label-free: does LOCO TPM renormalisation (exclude chr(A) transcripts) remove the negative-t skew?
Recomputes M1 t for G1-passing regulators on chr 1 and 13 with the original and renormalised Y."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/Users/colmclaude/ClaudeExperiments/DepMapTR/src')
from depmaptr import io as dio
from depmaptr.pipeline import Params, prepare, GroupCache, group_t
from depmaptr.assoc import centred_norm
ds = dio.load_real('S2', checksums=False)
p = Params()
prep = prepare(ds, p)
z = np.load('/Users/colmclaude/ClaudeExperiments/DepMapTR/data/results/scores/g1_t100_w0__a3_r095_s2_w0/scores.npz')
regs = z['regulators'].astype(str); cis = dict(zip(regs, z['cis_t']))
echrom = ds.genes.loc[ds.expr_keys, 'chrom'].astype(str).to_numpy()
tchrom = ds.genes.loc[ds.targets, 'chrom'].astype(str).to_numpy()
tpm = 2.0 ** ds.expr - 1.0
tot = tpm.sum(1)
YT0 = prep.YT.copy(); te = np.array([prep.ex_idx[k] for k in ds.targets])
mu = ds.expr[:, te].mean(0); q = pd.qcut(mu, 5, labels=False)
rng = np.random.default_rng(0)
for c in ['1', '13', '8']:
    keys = [k for k in ds.regulators if ds.genes.loc[k, 'chrom'] == c and cis[k] > 0 and cis[k]**2 >= 100]
    keys = np.array(rng.choice(keys, min(150, len(keys)), replace=False))
    share = tpm[:, echrom == c].sum(1) / tot
    out = {}
    for nm in ['orig', 'renorm']:
        if nm == 'orig':
            prep.YT = YT0
        else:
            prep.YT = np.log2(tpm[:, te] / (1.0 - share)[:, None] + 1.0)
        prep.y_raw = centred_norm(prep.YT)
        cache = GroupCache(prep)
        gv = cache.get(c)
        t, df, _ = group_t(prep, gv, keys)
        tt = t[:, tchrom != c]
        out[nm] = tt
        rm = np.nanmean(tt, 1)
        qs = ' '.join(f'{np.nanmean(tt[:, (q[tchrom != c] == k)]):+.3f}' for k in range(5))
        print(f'chr{c} {nm:6s}: n regs {len(keys)}, row-mean t median {np.median(rm):+.3f}, frac t<=-3 {np.nanmean(tt <= -3):.4f}, frac t>=3 {np.nanmean(tt >= 3):.4f}; mean t by target-expr quintile: {qs}')
    r = [np.corrcoef(np.nan_to_num(out['orig'][i]), np.nan_to_num(out['renorm'][i]))[0, 1] for i in range(len(keys))]
    print(f'   corr(orig row, renorm row) median {np.median(r):.3f}')
