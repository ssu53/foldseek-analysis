# %%

import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt



# Read the raw tmalign output (compact tabular format)
# skip top 5 rows where you left annotations

tmalign_raw = pd.read_csv(
    # 'outputs_2024-11-30-21-55-05/tmalign.out', # family
    # 'outputs_2024-11-30-22-14-44/tmalign.out', # superfamily
    # 'outputs_2024-11-30-22-26-34/tmalign.out', # fold

    # 'outputs_2024-12-01-00-08-03/tmalign.out',        # done, within fold, 56418535
    # 'outputs_2024-12-01-00-08-12/tmalign.out',        # done, split_00, 56419044
    # 'outputs_2024-12-01-00-51-24-26600/tmalign.out',  # done, split_05, 56420142
    # 'outputs_2024-12-01-01-13-58-70832/tmalign.out',  # done, split_02, 56421842
    # 'outputs_2024-12-01-01-25-14-40274/tmalign.out',  # done, split_04, 56422211
    # 'outputs_2024-12-01-01-55-55-75521/tmalign.out',  # done, split_01, 56422780
    # 'outputs_2024-12-01-02-09-30-72950/tmalign.out',  # done, split_03, 56423078

    # 'outputs_2024-12-05-14-23-40-81341/tmalign.out', # val 4k
    # 'outputs_2024-12-05-14-26-57-50236/tmalign.out', # val 2k above 0.6
    # 'outputs_2024-12-03-18-18-02-60350/tmalign.out', # val 2k random pairs
    # 'outputs_2024-12-03-18-22-00-55645/tmalign.out', # val 2k within fold pairs
    'outputs_2024-12-03-18-38-28-78170/tmalign.out', # val all within fold

    # 'outputs_2024-12-05-13-22-17-64530/tmalign.out', # train 20k above tms 0.6
    
    
    sep='\t',
    header=None,
    skiprows=5,
    # skipfooter=1, # this hack makes it use python engine
    engine='python',
    )




cpu_time = tmalign_raw.loc[2::3,:].dropna(axis=1)[0]
cpu_time = cpu_time.apply(lambda x: x.replace('Total CPU time is ','').replace(' seconds',''))
cpu_time = cpu_time.astype(float)
cpu_time = cpu_time.reset_index(drop=True)



tmalign = tmalign_raw.loc[1::3,:].dropna(axis=1)
tmalign.columns = [
    'prot_1',
    'prot_2',
    'tms_1',
    'tms_2',
    'rmsd',
    'seq_id_1',
    'seq_id_2',
    'seq_id',
    'len_1',
    'len_2',
    'len_aln',
]
tmalign = tmalign.apply(pd.to_numeric, errors='ignore')

path = '/scratch/groups/jamesz/shiye/scope40/'
tmalign.prot_1 = tmalign.prot_1.apply(lambda x: x.replace(path,''))
tmalign.prot_2 = tmalign.prot_2.apply(lambda x: x.replace(path,''))
tmalign.reset_index(drop=True, inplace=True)

tmalign['cpu_time'] = cpu_time
tmalign['tms'] = tmalign[['tms_1','tms_2']].mean(axis=1).round(4)

print(f"num samples {len(tmalign)}")
print(f"mean cpu time {tmalign.cpu_time.sum() / len(tmalign) :.5f}")


tmalign

# %%

tmalign[['prot_1','prot_2','tms']].to_csv(
    '/home/groups/jamesz/shiye/3d_protein_probing/data/embed_for_retrieval/train_data/tmaln_data_val_2k.csv',
    # '/home/groups/jamesz/shiye/protein_vector_retrieve/tmalign/tmaln.csv',
    sep=' ',
    header=None,
    index=None,
)

# %%

tmalign[['prot_1','prot_2','tms']].sample(20000, random_state=0).to_csv(
    '/home/groups/jamesz/shiye/3d_protein_probing/data/embed_for_retrieval/train_data/tmaln_data_train_within_fold_20k.csv',
    sep=' ',
    header=None,
    index=None,
)

# %%

tmalign[['prot_1','prot_2','tms']].sample(2000, random_state=0).to_csv(
    '/home/groups/jamesz/shiye/3d_protein_probing/data/embed_for_retrieval/train_data/tmaln_data_val_within_fold_2k.csv',
    sep=' ',
    header=None,
    index=None,
)
# %%

plt.figure()
plt.hist(tmalign.tms, bins=20)
plt.show()

# %%

tmalign = tmalign_raw[(~tmalign_raw[0].str.startswith('Total CPU time')) & (tmalign_raw[0] != '#PDBchain1')]
tmalign


# %%
