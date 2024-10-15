# %%

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# %%

# TODO: create alignment files based on TMscores, etc. to eval them on fold classif.

df = pd.read_csv('../training/tmp/result.rocx', sep='\t')

for level in ['FAM', 'SFAM', 'FOLD']:

    sensitivities = np.linspace(0., 1., 20)
    cumul = []
    for sens in sensitivities:
        cumul.append((df[level] >= sens).mean().item())

    plt.figure(figsize=(5,3))
    plt.scatter(cumul, sensitivities)
    plt.grid()
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.title(level)
    plt.xlabel('Fraction of queries')
    plt.ylabel('Sensitivity up to the 1st FP')
    plt.show()

# %%
