# %%

from itertools import combinations
import numpy as np
import pandas as pd

# %%


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



def main():

    split = 'train'
    num_pairs = 20000

    valid_pdb_ids = pd.read_csv(
        f'/home/groups/jamesz/shiye/foldseek-analysis/training/data_dev/valid_pdb_ids_{split}.csv',
        header=None,
        )
    valid_pdb_ids = valid_pdb_ids[0]

    get_random_pairs(valid_pdb_ids, split, num_pairs=num_pairs)
    # get_random_subset(valid_pdb_ids, split)



if __name__ == '__main__':
    main()
