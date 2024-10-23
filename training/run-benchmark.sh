#!/bin/bash -x

ENCODER=$1
STATES=$2
SUBMAT=$3

PDBS=$4 # list of SIDs
SCOPLOOKUP=$5

THETA=$6 # virtc
TAU=$7
D=$8

INVALIDSTATES=$9

MODE=${10}

PDB_DIR='/oak/stanford/groups/jamesz/shiye/scope40/'

# Filter scop_lookup.tsv
awk 'FNR==NR {pdbs[$1]=1; next}
     ($1 in pdbs) {print $0}' \
        $PDBS $SCOPLOOKUP > tmp/scop_lookup_filtered.tsv

# Encode validation PDBs
$RUN \
  ./encode_pdbs.py $ENCODER $STATES \
  --pdb_dir $PDB_DIR --virt $THETA $TAU $D \
  --invalid-state $INVALIDSTATES \
  --encoder_feat_type $MODE \
  < $PDBS > tmp/seqs_.csv

cp $SUBMAT tmp/sub_score.mat

# Prepare for benchmark
mkdir -p tmp/splits
mkdir -p tmp/alignments
awk '{print ">" $1} {print $2}' < tmp/seqs_.csv > tmp/target.fasta
split -n 30 -d tmp/target.fasta tmp/splits/split_ --additional-suffix=.fasta

./run-smithwaterman.sh 8 2  # sub_score.mat, target.fasta

# buffering issue (?) wait for all the files produced from the above script to appear
sleep 30 
# wait

./roc1.awk tmp/scop_lookup_filtered.tsv \
    <(cat tmp/alignments/*.m8) > tmp/result.rocx

# Calculate AUC
awk '{famsum+=$3; supfamsum+=$4; foldsum+=$5} END{print famsum/NR,supfamsum/NR,foldsum/NR}' \
    tmp/result.rocx
