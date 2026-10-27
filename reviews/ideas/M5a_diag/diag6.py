"""Label-free: global sign skew of trans t vs target expression level (compositional TPM artefact?)."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
R = '/Users/colmclaude/ClaudeExperiments/DepMapTR/data/results/scores/'
z = np.load(R + 'g1_t100_w0__a3_r095_s2_w0/scores.npz')
regs = z['regulators'].astype(str); tg = z['targets'].astype(str); cis = z['cis_t']
ann = pd.read_csv('/Users/colmclaude/ClaudeExperiments/DepMapTR/data/processed/gene_annotation.tsv', sep='\t', dtype={'chrom': str}).set_index('depmap_column')
rc = ann.loc[regs, 'chrom'].values; tc = ann.loc[tg, 'chrom'].values
g1 = (cis > 0) & (cis**2 >= 100)
T = z['t'][g1]; rcg = rc[g1]
mask = rcg[:, None] != tc[None, :]
Tm = np.where(mask, T, np.nan)
rowmean = np.nanmean(Tm, 1)
print('G1 regs n=%d: row-mean trans t: median %.3f, q10 %.3f, q90 %.3f, frac<0 %.2f' % (len(rowmean), np.median(rowmean), *np.percentile(rowmean, [10, 90]), (rowmean < 0).mean()))
print('global: frac t>=3 %.4f, frac t<=-3 %.4f' % (np.nanmean(Tm >= 3), np.nanmean(Tm <= -3)))
E = pd.read_parquet('/Users/colmclaude/ClaudeExperiments/DepMapTR/data/processed/depmap_expr_matched.parquet').set_index('ModelID')
panel = pd.read_csv('/Users/colmclaude/ClaudeExperiments/DepMapTR/data/processed/panel_s2.tsv', sep='\t')
idc = [c for c in panel.columns if c.lower() in ('modelid', 'model_id')][0]
E = E.loc[panel[idc]]
mu = E[tg].mean(0).values; sd = E[tg].std(0).values
colmean = np.nanmean(Tm, 0)
print('target column-mean t (over G1 regs) vs target mean expr S2: spearman %.3f; vs SD %.3f; vs mean/SD %.3f' % (spearmanr(colmean, mu)[0], spearmanr(colmean, sd)[0], spearmanr(colmean, mu / np.maximum(sd, 1e-3))[0]))
q = pd.qcut(mu, 5, labels=False)
for k in range(5): print('  expr quintile %d: mean col-mean t %.3f, frac t<=-3 %.4f, frac t>=3 %.4f' % (k + 1, colmean[q == k].mean(), np.nanmean(Tm[:, q == k] <= -3), np.nanmean(Tm[:, q == k] >= 3)))
# compositional covariate: TPM share of each chromosome per line
tpm = 2 ** E.values - 1
allc = ann.loc[E.columns, 'chrom'].values
tot = tpm.sum(1)
by = {}
for c in ['1', '8', '20', '13', '12']:
    share = tpm[:, allc == c].sum(1) / tot
    by[c] = share
    # regulators on c: corr of row-mean t with corr(x_A, share)? Simple: corr of arm-level share with total
    print(f'chr{c}: TPM share median {np.median(share):.3f}, CV {share.std()/share.mean():.3f}')
# relate per-regulator row mean to corr(x_A, log share of chr(A)) using CN
CN = pd.read_parquet('/Users/colmclaude/ClaudeExperiments/DepMapTR/data/processed/depmap_cn_matched.parquet').set_index('ModelID').loc[panel[idc]]
rg = regs[g1]
sel = [i for i, r in enumerate(rg) if rcg[i] in by]
cc = []
for i in sel:
    x = np.log2(np.clip(CN[rg[i]].values, 0.25, 8)); s = np.log(by[rcg[i]])
    ok = np.isfinite(x) & np.isfinite(s); cc.append(np.corrcoef(x[ok], s[ok])[0, 1])
print('G1 regs on chr 1/8/20/13/12 (n=%d): spearman(row-mean t, corr(x_A, log TPM-share of chr(A))) = %.3f; median corr(x_A, share) %.3f' % (len(sel), spearmanr(rowmean[sel], cc)[0], np.median(cc)))
