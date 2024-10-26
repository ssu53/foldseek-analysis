# %%

import sys
import os
import numpy as np
import pandas as pd



def get_alignment_files(
    alignments_dir = '/home/groups/jamesz/shiye/foldseek-analysis/training/tmp19_v4_encodings/alignments'
):
    alignment_files = os.listdir(alignments_dir)
    alignment_files = sorted(alignment_files)
    alignment_files = [os.path.join(alignments_dir, fn) for fn in alignment_files]
    return alignment_files



def calculate_num_samples_per_file(alignment_files, total_samples):

    num_pairs_per_file = np.full(len(alignment_files), None)

    # Get the number of pairs in each file

    for i,alignment_file in enumerate(alignment_files):

        with open(alignment_file) as f:
            num_pairs_per_file[i] = sum(1 for line in f if line.rstrip('\n'))

    print(f"{num_pairs_per_file=}")

    # Sample that many from each

    num_samples_per_file = np.ceil(num_pairs_per_file / np.sum(num_pairs_per_file) * total_samples)
    print(f"{num_samples_per_file=}")
    print(f"{np.sum(num_samples_per_file)=}")

    return num_samples_per_file



def get_samples_per_file(alignment_files, num_samples_per_file, seed=0):
        
    pairs_all = pd.DataFrame()

    for i,alignment_file in enumerate(alignment_files):

        # Read the pair alignment file
        pairs = pd.read_csv(alignment_file, header=None, sep=' ',)
        
        # Sample from it
        pairs = pairs.sample(num_samples_per_file[i], random_state=seed)

        # Update the global dataframe of pairs
        pairs_all = pd.concat((pairs_all, pairs), ignore_index=True)
    

    pairs_all.index = pairs_all.apply(lambda x: f"{x[0]}-{x[1]}", axis=1)
    assert len(pairs_all) == np.sum(num_samples_per_file)
    # print(pairs_all)
    # print(num_samples_per_file)
    # print(np.sum(num_samples_per_file))

    return pairs_all



def filter_already_computed_samples(
    pairs_all,
    master_file = 'tmalign.csv',
):

    # Filter pairs_all for those which are not already calculated 

    already_computed_pairs = pd.read_csv(master_file, sep='\t', usecols=['prot_1','prot_2'])
    already_computed_pairs = already_computed_pairs.apply(lambda x: f"{x.prot_1}-{x.prot_2}", axis=1)
    already_computed_pairs = set(already_computed_pairs.tolist())

    pairs_all_filtered = pairs_all[~pairs_all.index.isin(already_computed_pairs)]
    print(f"{len(pairs_all_filtered)=}")
    
    return pairs_all_filtered



def main(alignments_dir, num_samples, master_file, pair_file):

    alignment_files = get_alignment_files(alignments_dir)
    num_samples_per_file = calculate_num_samples_per_file(alignment_files, num_samples)
    pairs = get_samples_per_file(alignment_files, num_samples_per_file)
    pairs = filter_already_computed_samples(pairs, master_file)
    pairs.to_csv(pair_file, index=None, header=None, sep=' ')


if __name__ == '__main__':
    alignments_dir = sys.argv[1]
    num_samples = int(sys.argv[2])
    master_file = sys.argv[3]
    pair_file = sys.argv[4]
    main(alignments_dir, num_samples, master_file, pair_file)
