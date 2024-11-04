# %%

import os
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt

# %%

alignments_path = '../training/tmp0_foldseek_filtered/alignments'
# alignments_path = '../training/tmp19_v4_encodings/alignments'
# alignments_path = '../training/tmp21_v6_encodings/alignments'

metrics_path = 'MASTER_RESULTS_v0.csv'

metrics = pd.read_csv(metrics_path, sep='\t')
sid_pair_to_metrics_ind = {f"{row.prot_1}-{row.prot_2}": ind for ind,row in metrics.iterrows()}


alignment_fns = os.listdir(alignments_path)
alignment_fns = sorted(alignment_fns)

# %%


found_count = 0

for alignment_fn in alignment_fns:

    print(alignment_fn)

    alignments = pd.read_csv(os.path.join(alignments_path, alignment_fn), header=None, sep=' ')
    alignments.columns = ['sid1', 'sid2', 'score']
    # alignments = {f"{row.sid1}-{row.sid2}": row.score for index,row in alignments.iterrows()}

    for i in tqdm(alignments.index):

        key = f"{alignments.sid1[i]}-{alignments.sid2[i]}"
        if key in sid_pair_to_metrics_ind:
            found_count += 1
            j = sid_pair_to_metrics_ind[key]
            metrics.loc[j, 'new'] = alignments.score[i]

    print(f"So far found {found_count}.")

# %%



def plot_and_correlate(metrics, xvar, yvar, corr_method='spearman', top_only=False, clean_outliers=False, verbose=False):

    if top_only:
        if xvar in {'tms_1', 'tms_2'}:
            metrics = metrics[metrics[xvar] > 0.6]
        elif xvar == 'new':
            metrics = metrics[metrics[xvar] > 200]
        else:
            raise NotImplementedError
    
    if clean_outliers:
        outlier_pctl = 98
        if xvar in {'rmsd', 'lddt', 'chamfer', 'emd'}:
            x_cutoff = metrics[xvar].describe(percentiles=[outlier_pctl/100])[f'{outlier_pctl}%']
            metrics = metrics[metrics[xvar] < x_cutoff]
        if yvar in {'rmsd', 'lddt', 'chamfer', 'emd'}:
            y_cutoff = metrics[yvar].describe(percentiles=[outlier_pctl/100])[f'{outlier_pctl}%']
            metrics = metrics[metrics[yvar] < y_cutoff]
        # if xvar in {'tms_1', 'tms_2', 'new'}:
        #     raise NotImplementedError # trim from bottom
        # if yvar in {'tms_1', 'tms_2', 'new'}:
        #     raise NotImplementedError # trim from bottom

    if verbose:
        fig = plt.figure(figsize=(4,4))
        plt.scatter(metrics[xvar], metrics[yvar], s=5)
        plt.xlabel(xvar)
        plt.ylabel(yvar)
        plt.show()

    corr = metrics[xvar].corr(metrics[yvar], method=corr_method).item()
    # print(f"correlation: {corr:.5f}")

    return corr

# %%

print(f"{len(metrics)} rows")
metrics = metrics[~metrics.new.isna()]
print(f"{len(metrics)} data points")
metrics = metrics[metrics.rmsd > 0.01] # remove those that are like 0 which make correlation misleading
print(f"{len(metrics)} data points (after removing near identical matches...)")

# %%

df_corrs_to_tms = pd.DataFrame()
df_corrs_to_tms.index.name = 'tms_2'
for i,metric in enumerate(['rmsd', 'lddt', 'chamfer', 'emd', 'new']):
    df_corrs_to_tms.loc[metric, 'all'] = plot_and_correlate(metrics, 'tms_2', metric, top_only=False, clean_outliers=metric=='chamfer', verbose=True)
    df_corrs_to_tms.loc[metric, 'top'] = plot_and_correlate(metrics, 'tms_2', metric, top_only=True, clean_outliers=metric=='chamfer')
display(df_corrs_to_tms)

# %%

df_corrs_to_new = pd.DataFrame()
df_corrs_to_new.index.name = 'new'
for metric in ['rmsd', 'lddt', 'chamfer', 'emd', 'tms_2']:
    df_corrs_to_new.loc[metric, 'all'] = plot_and_correlate(metrics, 'new', metric, top_only=False, clean_outliers=metric=='chamfer', verbose=True)
    df_corrs_to_new.loc[metric, 'top'] = plot_and_correlate(metrics, 'new', metric, top_only=True, clean_outliers=metric=='chamfer', verbose=True)
display(df_corrs_to_new)


# %%
