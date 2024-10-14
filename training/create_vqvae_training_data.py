#! /usr/bin/env python3
"""
    ./create_vqvae_training_data.py data/pdbs 270 0 2
"""

import os.path
import numpy as np
import random
import sys
from functools import partial

import extract_pdb_features
import util
from tqdm import tqdm

import torch
import esm

data_dir = os.path.join(os.path.dirname(__file__), 'data/')


feature_cache = {}  # path: (features, mask)

def encoder_features(pdb_path, virt_cb):
    """
    Calculate 3D descriptors for each residue of a PDB file.
    """
    feat = feature_cache.get(pdb_path, None)
    if feat is not None:
        return feat

    coords, valid_mask = extract_pdb_features.get_coords_from_pdb(pdb_path, full_backbone=True)
    coords = extract_pdb_features.move_CB(coords, virt_cb=virt_cb)
    partner_idx = extract_pdb_features.find_nearest_residues(coords, valid_mask)
    features, valid_mask2 = extract_pdb_features.calc_angles_forloop(coords, partner_idx, valid_mask)

    seq_dist = (partner_idx - np.arange(len(partner_idx)))[:, np.newaxis]
    log_dist = np.sign(seq_dist) * np.log(np.abs(seq_dist) + 1)

    vae_features = np.hstack([features, log_dist])
    feature_cache[pdb_path] = vae_features, valid_mask2

    return vae_features, valid_mask2


def align_features(pdb_dir, virtual_center, sid1, sid2, cigar_string):
    """
    Return aligned descriptors for a given alignment between two PDBs.
    """

    idx_1, idx_2 = util.parse_cigar(cigar_string).T

    feat1, mask1 = encoder_features(os.path.join(pdb_dir, sid1), virtual_center)
    feat2, mask2 = encoder_features(os.path.join(pdb_dir, sid2), virtual_center)

    valid_mask = mask1[idx_1] & mask2[idx_2]
    idx_1 = idx_1[valid_mask]
    idx_2 = idx_2[valid_mask]

    x = np.vstack([feat1[idx_1], feat2[idx_2]])
    y = np.vstack([feat2[idx_2], feat1[idx_1]])
    return x, y  # (n x 10, n x 10)


def encoder_features_position(pdb_path, virt_cb):
    """
    Return the virtual center positions of each residue
    """
    feat = feature_cache.get(pdb_path, None)
    if feat is not None:
        return feat
    
    # Read the coordinates from PDB file. Returns coordinates of CA + CB + N + C.
    coords, valid_mask = extract_pdb_features.get_coords_from_pdb(pdb_path, full_backbone=True)
    assert coords.shape[1] == 4 * 3 

    # Replace CB coordinate with coordinate of virtual center. Note that move_CB changes coords in place!
    coords = extract_pdb_features.move_CB(coords, virt_cb=virt_cb)
    coords_cb = coords[:, 3:6]

    features_cache[pdb_path] = coords_cb, valid_mask

    return coords_cb, valid_mask


def encoder_features_esm(pdb_path, model, alphabet, batch_converter, device):

    feat = feature_cache.get(pdb_path, None)
    if feat is not None:
        return feat

    coords, valid_mask = extract_pdb_features.get_coords_from_pdb(pdb_path, full_backbone=True)

    prot_name = pdb_path.split('/')[-1]

    with open(f"tmp/fasta/{prot_name}.fasta", 'r') as f:
        fasta_name, fasta_seq = f.readlines()
    fasta_name = fasta_name.strip().replace('>','')
    fasta_seq = fasta_seq.strip()
    data = [(fasta_name, fasta_seq)]
    assert prot_name == fasta_name
    assert coords.shape[0] == len(fasta_seq)

    # Prepare data
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
    batch_tokens = batch_tokens.to(device)

    # Residue level representations
    FINAL_LAYER = 6
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[FINAL_LAYER], return_contacts=True)
    token_representations = results["representations"][FINAL_LAYER]
    token_representations = token_representations.squeeze(0)
    token_representations = token_representations[1 : batch_lens - 1]

    # print(batch_lens, len(fasta_seq), token_representations.shape, valid_mask.shape)
    assert len(fasta_seq) == len(valid_mask)
    assert token_representations.shape == (len(fasta_seq), 320)
    token_representations = token_representations.cpu().numpy()

    features_cache[pdb_path] = token_representations, valid_mask

    return token_representations, valid_mask


def align_features_esm(pdb_dir, sid1, sid2, cigar_string):
    """
    Return aligned ESM descriptors for a given alignment between two PDBs.
    """

    idx_1, idx_2 = util.parse_cigar(cigar_string).T

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    model.to(device)

    feat1, mask1 = encoder_features_esm(os.path.join(pdb_dir, sid1), model, alphabet, batch_converter, device)
    feat2, mask2 = encoder_features_esm(os.path.join(pdb_dir, sid2), model, alphabet, batch_converter, device)
    # print(f"{idx_1.shape=} {idx_2.shape=} {feat1.shape=} {mask1.shape=} {feat2.shape=} {mask2.shape=}")

    valid_mask = mask1[idx_1] & mask2[idx_2]
    idx_1 = idx_1[valid_mask]
    idx_2 = idx_2[valid_mask]

    x = np.vstack([feat1[idx_1], feat2[idx_2]])
    y = np.vstack([feat2[idx_2], feat1[idx_1]])
    return x, y 


def encoder_features_custom(pdb_id, pdb_id_to_encoding, pdb_id_to_protein):
    """
    """

    # Get encoding
    encoding = pdb_id_to_encoding[pdb_id]
    assert encoding.ndim == 2
    assert encoding.size(1) == 10 # specific to the current encoding

    # Get valid_mask from protein object
    mask = pdb_id_to_protein[pdb_id]['valid_mask']
    assert mask.shape == (encoding.size(0),)

    # Mark the mask positions at end residues as False
    # This is consistent with extract_pdb_features.calc_angles_forloop
    # due to incomplete angular neighbourhood descriptors
    mask[0] = False
    mask[-1] = False

    return encoding.numpy(), mask.numpy()


def align_features_custom(pdb_dir, sid1, sid2, cigar_string, pdb_id_to_encoding, pdb_id_to_protein):
    """
    Return aligned custom descriptors for a given alignment between two PDBs.
    """

    idx_1, idx_2 = util.parse_cigar(cigar_string).T

    feat1, mask1 = encoder_features_custom(sid1, pdb_id_to_encoding, pdb_id_to_protein)
    feat2, mask2 = encoder_features_custom(sid2, pdb_id_to_encoding, pdb_id_to_protein)

    # # check they are the same shape as your entries in proteins.pt
    # pdb_id_to_protein = torch.load('../../3d_protein_probing/data/scope40_foldseek_compatible/proteins_train.pt')
    # structure1 = pdb_id_to_protein[sid1]['structure']
    # structure2 = pdb_id_to_protein[sid2]['structure']
    # coords1, valid_mask1 = extract_pdb_features.get_coords_from_pdb(os.path.join(pdb_dir, sid1), full_backbone=True)
    # coords2, valid_mask2 = extract_pdb_features.get_coords_from_pdb(os.path.join(pdb_dir, sid2), full_backbone=True)
    # assert len(coords1) == len(feat1) == len(valid_mask1) == len(mask1)
    # assert len(coords2) == len(feat2) == len(valid_mask2) == len(mask2)

    valid_mask = mask1[idx_1] & mask2[idx_2]
    idx_1 = idx_1[valid_mask]
    idx_2 = idx_2[valid_mask]

    x = np.vstack([feat1[idx_1], feat2[idx_2]])
    y = np.vstack([feat2[idx_2], feat1[idx_1]])
    return x, y  # (n x 10, n x 10)


if __name__ == '__main__':
    pdb_dir = sys.argv[1]
    pairfile = sys.argv[2]
    a = int(sys.argv[3])
    b = int(sys.argv[4])
    c = float(sys.argv[5])
    out = sys.argv[6]
    mode = sys.argv[7]
    virtual_center = (a, b, c)

    with open(data_dir + 'pdbs_train.txt') as file:
        pdbs_train = set(file.read().splitlines())

    # Find alignments between PDBs of the training set
    alignments = []
    with open(pairfile) as file:
        for line in file:
            sid1, sid2, cigar_string = line.rstrip('\n').split()

            if sid1 in pdbs_train and sid2 in pdbs_train:
                #print(' '.join((sid1, sid2, cigar_string)))
                alignments.append((sid1, sid2, cigar_string))

    # Not needed execept to exactly reproduce result
    random.Random(123).shuffle(alignments)

    print(f"{len(alignments)=}")

    if mode == 'foldseek-default':
        align_features_fn = partial(align_features, pdb_dir=pdb_dir, virtual_center=virtual_center)
    elif mode == 'esm':
        align_features_fn = partial(align_features_esm, pdb_dir=pdb_dir)
    elif mode == 'custom':
        pdb_id_to_encoding = torch.load('data_dev/encodings.pt')
        pdb_id_to_protein = torch.load('data_dev/proteins.pt')
        align_features_fn = partial(
            align_features_custom,
            pdb_dir=pdb_dir,
            pdb_id_to_encoding=pdb_id_to_encoding,
            pdb_id_to_protein=pdb_id_to_protein,
        )
    else:
        raise NotImplementedError
    
    xy = []  # (n x 10, n x 10)
    for sid1, sid2, cigar_string in tqdm(alignments):
        xy.append(align_features_fn(sid1=sid1, sid2=sid2, cigar_string=cigar_string))

    # Write features to disc
    x_feat = np.vstack([x for x, y in xy])
    y_feat = np.vstack([y for x, y in xy])
    idx = np.arange(len(x_feat))
    np.random.RandomState(123).shuffle(idx)

    print(f"{x_feat.shape=} {y_feat.shape=}")

    np.save(out, np.dstack([x_feat[idx], y_feat[idx]]))

