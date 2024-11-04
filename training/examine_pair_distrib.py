# %%

from tqdm import tqdm

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import (
    mean_absolute_error,
    mean_absolute_percentage_error,
    mean_squared_error,
    r2_score,
)
from sklearn.decomposition import PCA

# %%

neural_v0 = np.load('training/tmp3_v0_encodings/vaevq_training_data.npy')
neural_v1 = np.load('training/tmp4_v1_encodings/vaevq_training_data.npy')
foldseek_v1 = np.load('training/tmp5_foldseek_on_valid_ids/vaevq_training_data.npy')
# foldseek_v1_ = np.load('training/tmp6_foldseek_ignorehetatm/vaevq_training_data.npy') # this is allclose to above

# %%

print(f"{neural_v0.shape=}")
print(f"{neural_v1.shape=}")
print(f"{foldseek_v1.shape=}")

# %%

vaevq_training_data = {
    'foldseek_v1': foldseek_v1,
    'neural_v0': neural_v0,
    'neural_v1': neural_v1,
}
summary_stats = {}

# %%

num_samples_to_plot = 1000
verbose = False
do_pca = False

for key in vaevq_training_data:

    if key not in {
        'foldseek_v1',
        # 'neural_v0',
        'neural_v1',
    }: 
        continue
    
    print(key)
    
    data = vaevq_training_data[key]
    num_feats = data.shape[1]
    stats = pd.DataFrame(
        index=range(num_feats),
        columns=['mean','std','rmse','r2'])
    summary_stats[f"{key}{'_pca' if do_pca else ''}"] = stats

    data_x = data[:,:,0]
    data_y = data[:,:,1]

    if do_pca:
        pca = PCA(n_components=num_feats)
        pca.fit(data_x)
        data_x = pca.transform(data_x)
        data_y = pca.transform(data_y)


    for i_feat in range(num_feats):

        feat_x = data_x[:, i_feat]
        feat_y = data_y[:, i_feat]
        assert feat_x.shape == feat_y.shape == (len(data),)
        # rmse, mae, mape, r2, plot, std

        feat_mean = np.mean(feat_x)
        feat_std = np.std(feat_x)

        assert np.allclose(feat_mean, np.mean(feat_y))
        assert np.allclose(feat_std, np.std(feat_y))

        feat_rmse = mean_squared_error(feat_x, feat_y)
        feat_r2 = r2_score(feat_x, feat_y)

        stats.loc[i_feat, 'mean'] = feat_mean
        stats.loc[i_feat, 'std'] = feat_std
        stats.loc[i_feat, 'rmse'] = feat_rmse
        stats.loc[i_feat, 'r2'] = feat_r2

        if verbose:

            print(f"feature {i_feat}: {feat_mean:.3f} +/- {feat_std:.3f}")
            print(f"RMSE {feat_rmse:.3f}")
            print(f"R2 {feat_r2:.3f}")

            ind_sampled = np.random.choice(range(len(data)), num_samples_to_plot)

            plt.figure()
            plt.scatter(feat_x[ind_sampled], feat_y[ind_sampled])
            plt.show()
    

# %%

summary_stats['foldseek_v1'].style.format('{:.3f}')

# %%

summary_stats['neural_v0'].style.format('{:.3f}')

# %%

summary_stats['neural_v1'].style.format('{:.3f}')

# %%

summary_stats['foldseek_v1_pca'].style.format('{:.3f}')

# %%
