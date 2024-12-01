# %%

from tqdm import tqdm
from itertools import combinations
import itertools
import random
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser



def get_random_pairs(
    valid_pdb_ids,
    split,
    num_pairs = 5000,
    seed = 0,
):
    """
    Random pairs sampling strategy.
    Not enough high-aligned hits.
    """

    np.random.seed(seed)

    pair_pdb_ids = np.random.choice(valid_pdb_ids, size=(num_pairs, 2))
    non_repeats = pair_pdb_ids[:,0] != pair_pdb_ids[:,1]
    pair_pdb_ids = pair_pdb_ids[non_repeats]
    print(f"{pair_pdb_ids.shape=}")


    pd.DataFrame(pair_pdb_ids).to_csv(
        f'/home/groups/jamesz/shiye/foldseek-analysis/metrics/pairfile_random_pairs_{split}.out',
        index=False,
        header=None,
        sep=' ',
        )



def get_random_subset(
    valid_pdb_ids,
    split,
    num_singles = 128,
    seed = 0,
):
    """
    Random subset, then all-vs-all (unique pair combinations).
    """
    np.random.seed(seed)

    pdb_id_samples = valid_pdb_ids.sample(num_singles)
    pair_pdb_ids = list(combinations(pdb_id_samples, r=2))
    print(f"{len(pair_pdb_ids)=}")

    pd.DataFrame(pair_pdb_ids).to_csv(
        f'/home/groups/jamesz/shiye/foldseek-analysis/metrics/pairfile_random_subset_{split}.out',
        index=False,
        header=None,
        sep=' ',
    )



def get_scop_lookup():
    
    scop_lookup = pd.read_csv(
        '/home/groups/jamesz/shiye/foldseek-analysis/training/data/scop_lookup.tsv',
        sep='\t', header=None)
    scop_lookup.columns = ['sid', 'family']
    scop_lookup.set_index('sid', inplace=True)
    scop_lookup['superfamily'] = scop_lookup['family'].apply(lambda x: x[:x.rfind('.')])
    scop_lookup['fold'] = scop_lookup['superfamily'].apply(lambda x: x[:x.rfind('.')])
    scop_lookup['class'] = scop_lookup['fold'].apply(lambda x: x[:x.rfind('.')])
    scop_lookup.sort_index(inplace=True)

    return scop_lookup



def check_pdb_num_models_and_num_chains():
    """
    found that for every pdb in valid_pdb_ids_train
    there is only one model and that model has only one chain
    """

    more_than_1_model = []
    more_than_1_chain = []

    for pdb_id in tqdm(valid_pdb_ids_train):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('None', f"/scratch/groups/jamesz/shiye/scope40/{pdb_id}")

        model = structure[0]  # take only first model
        chains = list(model.get_chains())
        # chain = list(model.get_chains())[0]  # take only first chain

        if len(structure) != 1:
            more_than_1_model.append(pdb_id)
            print(f"{pdb_id} {len(structure)} {len(chains)}")
        if len(chains) != 1:
            more_than_1_chain.append(pdb_id)
            print(f"{pdb_id} {len(structure)} {len(chains)}")
    
    print(f"{len(more_than_1_model)=}")
    print(f"{len(more_than_1_chain)=}")



def get_within_group_pairs(
    scop_lookup,
    valid_pdb_ids,
    level = 'fold',
    k = None,
    seed = 0,
):
    """
    :param level: family, superfamily, or fold
    :param k: number of samples in each group, if None takes all pairs
    :param seed: random seed for for samples
    """

    random.seed(seed)

    scop_lookup_subset = scop_lookup.loc[valid_pdb_ids,:]
    scop_lookup_subset.sort_index(inplace=True)
    scop_lookup_subset


    level_groups = scop_lookup_subset[level].unique()
    print(len(level_groups), level)


    all_pairs = []


    for group in tqdm(level_groups):

        pdb_ids_group = scop_lookup_subset.index[scop_lookup_subset[level] == group]
        pairs = list(itertools.combinations(pdb_ids_group, 2))

        if k is not None and len(pairs) > k:
            pairs = random.sample(pairs, k)

        all_pairs.extend(pairs)

    print(len(all_pairs))


    pd.DataFrame(all_pairs).to_csv(
        f'/home/groups/jamesz/shiye/foldseek-analysis/metrics/pairfile_within_{level}_train_{k if k is not None else "all"}.out',
        sep=' ', header=None, index=None)
    
    return all_pairs



def main1():

    split = 'train'
    num_pairs = 20000

    valid_pdb_ids = pd.read_csv(
        f'/home/groups/jamesz/shiye/foldseek-analysis/training/data_dev/valid_pdb_ids_{split}.csv',
        header=None,
        )
    valid_pdb_ids = valid_pdb_ids[0]

    get_random_pairs(valid_pdb_ids, split, num_pairs=num_pairs)
    # get_random_subset(valid_pdb_ids, split)



def main2():

    split = 'train'

    valid_pdb_ids = pd.read_csv(
        f'/home/groups/jamesz/shiye/foldseek-analysis/training/data_dev/valid_pdb_ids_{split}.csv',
        header=None,
    )[0].tolist()
    print(len(valid_pdb_ids))

    scop_lookup = get_scop_lookup()

    pairs_within_fold = get_within_group_pairs(scop_lookup, valid_pdb_ids, level='fold', k=None)


    # Get large set of pairs (~10% of all possible pairs, excluding those within fold)

    pairs = list(itertools.combinations(valid_pdb_ids, 2))
    pairs_sampled = random.sample(pairs, len(pairs) // 10)
    pairs_sampled_outside_fold = set(pairs_sampled) - set(pairs_within_fold)

    print(f"{len(pairs_within_fold)=}")
    print(f"{len(pairs)=}")
    print(f"{len(pairs_sampled)=}")
    print(f"{len(pairs_sampled_outside_fold)=}")

    pd.DataFrame(pairs_sampled_outside_fold).to_csv(
        f'/home/groups/jamesz/shiye/foldseek-analysis/metrics/pairfile_random_pairs_train_outside_fold.out',
        sep=' ', header=None, index=None)



if __name__ == '__main__':
    main2()

