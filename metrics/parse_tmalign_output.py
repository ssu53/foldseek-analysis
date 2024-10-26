# %%

import sys
import pandas as pd


def parse_tmalign_output(line, row, rows, pdb_dir):

    if line.strip() == '':
        return line, row, rows, pdb_dir

    elif line.startswith('### PDB_DIR'):
        _, _, pdb_dir = line.partition('PDB_DIR')
        pdb_dir = pdb_dir.strip()

    elif line.startswith(('###', '(You should use TM-score')):
        return line, row, rows, pdb_dir
    
    elif line.startswith("Name of Chain_1:"):
        _, _, after_keyword = line.partition("Name of Chain_1:")
        row.append(after_keyword.split()[0].replace(pdb_dir, ''))

    elif line.startswith("Name of Chain_2:"):
        _, _, after_keyword = line.partition("Name of Chain_2:")
        row.append(after_keyword.split()[0].replace(pdb_dir, ''))
    
    elif line.startswith("Length of Chain_1:"):
        _, _, after_keyword = line.partition("Length of Chain_1:")
        row.append(after_keyword.split()[0].replace(pdb_dir, ''))

    elif line.startswith("Length of Chain_2:"):
        _, _, after_keyword = line.partition("Length of Chain_2:")
        row.append(after_keyword.split()[0].replace(pdb_dir, ''))

    elif line.startswith("Aligned length"):
        line_split = line.split()
        assert len(line_split) == 7
        w0, w1, w2, w3, w4, w5, w6 = line_split
        assert w0 == 'Aligned'
        assert w1 == 'length='
        row.append(w2.replace(',',''))
        assert w3 == 'RMSD='
        row.append(w4.replace(',',''))
        assert w5 == "Seq_ID=n_identical/n_aligned="
        row.append(w6)
    
    elif line.startswith("TM-score="):
        _, _, after_keyword = line.partition("TM-score=")
        row.append(after_keyword.split()[0])
    
    elif line.startswith("(\":\" denotes residue pairs"):
        return line, row, rows, pdb_dir
    
    elif line.startswith("Total CPU time is"):
        _, _, after_keyword = line.partition("Total CPU time is")
        row.append(after_keyword.split()[0])
        # End of one pair alignment run output, pack into dict
        rows.append(row)
        row = []
    
    else:
        # The aligned residues
        row.append(line.replace('\n',''))

    return line, row, rows, pdb_dir



def get_parsed_df(load_path: str = 'tmalign.out') -> pd.DataFrame:

    pdb_dir = None
    rows = []
    row = []


    with open(load_path, 'r') as f:
        for line in f:
            line, row, rows, pdb_dir = parse_tmalign_output(line, row, rows, pdb_dir)
    print(f"{len(rows)=}")

    df = pd.DataFrame.from_dict(rows, orient='columns')
    df.columns = [
        'prot_1',
        'prot_2',
        'len_1',
        'len_2',
        'len_aln',
        'rmsd',
        'seq_id',
        'tms_1',
        'tms_2',
        'seq_1',
        'seq_aln',
        'seq_2',
        'cpu_time',
    ]
    df = df[[ # reorder
        'prot_1',
        'prot_2',
        'tms_1',
        'tms_2',
        'rmsd',
        'len_1',
        'len_2',
        'len_aln',
        'seq_id',
        'seq_1',
        'seq_aln',
        'seq_2',
        'cpu_time',
    ]]

    print(f"PDBs from {pdb_dir}")
    print(f"Ran and parsed {len(df)} pair alignments.")
    df.cpu_time = pd.to_numeric(df.cpu_time)
    print(f"Mean CPU time {df.cpu_time.mean():.3f}s ({df.cpu_time.sum():.3f}s total).")

    return df


def parse_seq_strings(seq1: str, seq2: str, seq_aln: str) -> str:
    cigar_str = []
    for s1, s2, sa in zip(seq1, seq2, seq_aln):
        if s1.isspace() and s2.isspace() and sa.isspace(): continue
        if sa == ':': cigar_str.append('P')
        elif sa == '.': cigar_str.append('M')
        elif sa == ' ' and s1 == '-': cigar_str.append('I')
        elif sa == ' ' and s2 == '-': cigar_str.append('D')
        else: 
            raise Exception(f"Invalid sequence strings {s1} {s2} {s1}")
    return ''.join(cigar_str)


def parse_cigar_string(string: str) -> str:
    if string == '': return string
    cigar_str = []
    char = string[0]
    cnt = 0
    for c in string:
        if c == char: 
            cnt += 1
        else: 
            cigar_str.append(str(cnt))
            cigar_str.append(char)
            cnt = 1
            char = c
    cigar_str.append(str(cnt))
    cigar_str.append(char)
    return ''.join(cigar_str)


def populate_cigar_strings(df: pd.DataFrame) -> pd.DataFrame:
    df['cigar_'] = df.apply(
        lambda row: parse_seq_strings(row.seq_1, row.seq_2, row.seq_aln),
        axis=1)
    df['cigar'] = df.apply(lambda row: parse_cigar_string(row.cigar_), axis=1)
    return df


def save_parsed_df_with_cigars(
    load_path: str = 'tmalign.out',
    save_path: str = 'tmalign.csv'
):

    df = get_parsed_df(load_path)
    df = populate_cigar_strings(df)

    if save_path is not None:
        df.to_csv(save_path, index=False, sep='\t')



def verify_against_file():
    """
    for pairs in tmaln-06-500.out, which contains TMalign outputs shipped with the repo
    """

    import numpy as np
    
    df1 = pd.read_csv('tmalign.csv', sep='\t')

    df2 = pd.read_csv('../training/data/tmaln-06_500.out', sep='\t', header=None)
    df2 = df2.dropna(how='all', axis=1)
    df2.columns = [
        'prot_1', 
        'prot_2',
        'tms',
        'tms_1',
        'tms_2',
        'rmsd',
        'len_1',
        'len_2',
        'len_aln',
        'cigar',
        ]

    assert all(df1.prot_1 == df2.prot_1)
    assert all(df1.prot_2 == df2.prot_2)
    assert np.allclose(df1.tms_1, df2.tms_1)
    assert np.allclose(df1.tms_2, df2.tms_2)
    assert all(df1.rmsd == df2.rmsd.round(2))
    assert all(df1.len_1 == df2.len_1)
    assert all(df1.len_2 == df2.len_2)
    assert all(df1.len_aln == df2.len_aln)
    assert all(df1.cigar == df2.cigar)

    print("Passed!")


if __name__ == '__main__':
    load_path = sys.argv[1] # raw outputs of TMalign
    save_path = sys.argv[2] # tabulated outputs with cigar strings
    save_parsed_df_with_cigars(load_path, save_path)
    # verify_against_file()

    

# %%
