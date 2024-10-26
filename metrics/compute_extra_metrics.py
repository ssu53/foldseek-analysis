# %% 

import sys
import re
from tqdm import tqdm
from pathlib import Path
from typing import List, Union

import pandas as pd
import numpy as np
import torch

from scipy.stats import wasserstein_distance_nd
from chamferdist import ChamferDistance
chamfer_dist_fn = ChamferDistance()

from lddt_dis import lddt, read_ca_from_pdb, cigar_to_tuples



def sym_chamfer_(
    point_cloud1: torch.Tensor,
    point_cloud2: torch.Tensor,
    distance_type: str = "vanilla",
) -> torch.Tensor:
    """
    Compute the symmetrised Chamfer Distance between two non-batched point clouds. No Nan values are allowed. Uses POT library
    :param point_cloud1: (P1, D)
    :param point_cloud2: (P2, D)
    :param distance_type: Either "vanilla" for the regular chamfer distance, or
        "density_aware" for the density-aware version (see https://arxiv.org/abs/2111.12702)
    :return: chamfer_dist (torch.Tensor): Shape (1) tensor representing the Chamfer Distance between the point clouds
    """
    # Compute pairwise distance matrices
    # This could also be done with cdist as below, but that's slower
    # dist_matrix = torch.cdist(point_cloud1, point_cloud2)
    dist_matrix = (point_cloud2[None, :, :] - point_cloud1[:, None, :]).pow(2).sum(-1)
    # Compute the minimum distance from points in point_cloud1 to point_cloud2 and vice versa
    forward_distances, indices1 = torch.min(dist_matrix, dim=-1)  # shape (P1,)

    backward_distances, indices2 = torch.min(dist_matrix, dim=-2)  # shape (P2,)

    if distance_type == "vanilla":
        chamfer_dist = 0.5 * (
            torch.mean(forward_distances) + torch.mean(backward_distances)
        )  # shape (1,)
    elif distance_type == "density_aware":
        exp_distances1 = torch.exp(-forward_distances)
        exp_distances2 = torch.exp(-backward_distances)

        # scale each distance by the number of indices that it matches
        weighted_distances1 = (
            exp_distances1
            / torch.bincount(indices1, minlength=point_cloud2.size(0)).float()[indices1]
        )
        weighted_distances2 = (
            exp_distances2
            / torch.bincount(indices2, minlength=point_cloud1.size(0)).float()[indices2]
        )

        chamfer_dist = 0.5 * (
            torch.mean(1 - weighted_distances1) + torch.mean(1 - weighted_distances2)
        )
    else:
        raise NotImplementedError(f"Distance type {distance_type} not implemented.")

    return chamfer_dist



def sym_chamfer(source_cloud, target_cloud):
    """
    Uses chamferdist package to compute Chamfer distances between two non-batched point clouds. 
    This runs faster than sym_chamfer_.
    Symmetrise by averaging forward and backwards.
    """
    assert source_cloud.ndim == target_cloud.ndim == 2
    dist_forward = chamfer_dist_fn(source_cloud[None,:,:], target_cloud[None,:,:])
    dist_backward = chamfer_dist_fn(target_cloud[None,:,:], source_cloud[None,:,:])
    dist = 0.5 * (dist_forward / len(source_cloud) + dist_backward / len(target_cloud))
    return dist



def emd(source_cloud, target_cloud):
    return wasserstein_distance_nd(source_cloud, target_cloud)



def read_tmalign_csv(
    path: Path,
    cols: List[str] = None,
    ) -> pd.DataFrame:
    
    df = pd.read_csv(path, sep='\t')

    if cols is None:
        df = df[[
            'prot_1',
            'prot_2',
            'tms_1',
            'tms_2',
            'rmsd',
            'len_1',
            'len_2',
            'len_aln',
            'cigar',
        ]]
    
    return df



def read_tmalign_transformation(
    dir: Path,
    fn: str,
):
    if isinstance(dir, str): dir = Path(dir)
    tf = np.genfromtxt(dir / fn, skip_header=2, skip_footer=7, usecols=(1,2,3,4))
    translation = tf[:,0]
    rotation = tf[:,1:4]

    return translation, rotation



def populate_metrics_for_row(
    df: pd.DataFrame,
    i: Union[str, int], 
    pdb_dir: Path,
    pdb_extension: str='',
    rot_mats_dir: str = 'rot_mats',
    compute_chamfer: bool = True,
    compute_emd: bool = True,
    compute_lddt: bool = True,
    ):
    
    # Get protein IDs
    id1 = df.prot_1[i]
    id2 = df.prot_2[i]
    cigar = df.cigar[i]

    # Get CA coordinates
    coords1 = read_ca_from_pdb(f"{pdb_dir}/{id1}{pdb_extension}")
    coords2 = read_ca_from_pdb(f"{pdb_dir}/{id2}{pdb_extension}")
    assert coords1.ndim == coords2.ndim == 3
    assert coords1.shape[0] == coords2.shape[0] == 1
    coords1 = coords1.squeeze(0)
    coords2 = coords2.squeeze(0)

    # Get sequence alignment
    aligned = cigar_to_tuples(cigar, query_start=0, target_start=0)
    
    # Get superposition coording to TMalign computed transformation
    translation, rotation = read_tmalign_transformation(dir=rot_mats_dir, fn=f"{id1}-{id2}.txt")
    coords1_tf = translation + (rotation @ coords1.T).T

    # LDDT
    # This metric is superposition-free; doesn't actually need the transformed coords1.
    # This implementation is batched. Unsqueeze.
    if compute_lddt:
        dist_lddt_per_residue = lddt(query_points=coords1_tf[None,:,:], target_points=coords2[None,:,:], aligned_pairs=aligned, per_residue=True)
        
        # Verify LDDT is superposition-free.
        # dist_lddt_per_residue_ = lddt(query_points=coords1[None,:,:], target_points=coords2[None,:,:], aligned_pairs=aligned, per_residue=True)
        # assert np.allclose(dist_lddt_per_residue, dist_lddt_per_residue_)

        df.loc[i,'lddt'] = np.round(dist_lddt_per_residue.mean(), 5)

    # Chamfer
    # This implementation requires input point clouds to be tensors.
    if compute_chamfer:
        with torch.no_grad():
            dist_chamfer = sym_chamfer(
                torch.tensor(coords1_tf).to(torch.float32),
                torch.tensor(coords2).to(torch.float32),
            ).item()
        df.loc[i,'chamfer'] = np.round(dist_chamfer, 5)

    # EMD
    if compute_emd:
        dist_emd = emd(coords1_tf, coords2)
        df.loc[i,'emd'] = np.round(dist_emd, 5)

    return df # has been modified in place



def main(
    load_path = 'tmalign.csv',
    pdb_dir = '/scratch/groups/jamesz/shiye/scope40',
    rot_mats_dir = 'rot_mats',
):

    df = read_tmalign_csv(load_path)

    for i in tqdm(df.index):
        populate_metrics_for_row(
            df,
            i,
            pdb_dir=pdb_dir,
            rot_mats_dir=rot_mats_dir,
            # compute_emd=False, # EMD computation is super slow.
        )

    return df


if __name__ == '__main__':
    load_path = sys.argv[1] # tabulated outputs of TMalign with cigar strings
    pdb_dir = sys.argv[2]   # directory of PDBs
    rot_mats_dir = sys.argv[3]  # directory of TMalign rotation matrices
    save_path = sys.argv[4] # tabulated outputs with extra metrics

    df = main(load_path, pdb_dir, rot_mats_dir)
    print(df[['rmsd','lddt','chamfer','emd']].mean())
    df = df[[ # reorder
        'prot_1',
        'prot_2',
        'tms_1',
        'tms_2',
        'rmsd',
        'len_1',
        'len_2',
        'len_aln',
        'lddt',
        'chamfer',
        'emd',
        'cigar',
    ]]
    df.to_csv(save_path, index=False, sep='\t')

