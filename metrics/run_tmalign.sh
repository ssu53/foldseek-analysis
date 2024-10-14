#!/bin/bash

# Set input arguments and defaults
PDB_DIR="${1:-$'/oak/stanford/groups/jamesz/shiye/scope40/'}"
PAIR_FILE="${2:-$'/home/groups/jamesz/shiye/foldseek-analysis/training/data/tmaln-06_500.out'}"
ROT_DIR="${3:-$'rot_mats/'}"
OUT_FILE="${4:-$'tmalign.out'}"


echo "Running..." $(date)

echo "Running TM-align (Version 20220412), compiled locally from cpp." > $OUT_FILE

echo "PDB_DIR" $PDB_DIR >> $OUT_FILE
echo "PAIR_FILE" $PAIR_FILE >> $OUT_FILE
echo "ROT_DIR" $ROT_DIR >> $OUT_FILE
echo "OUT_FILE" $OUT_FILE >> $OUT_FILE


# awk -v PDB_DIR="${PDB_DIR}" -v ROT_DIR="${ROT_DIR}" '{ system("echo " PDB_DIR $1 " " PDB_DIR $2 " -m " ROT_DIR $1 "-" $2 ".txt") }' $PAIR_FILE

# Run TM-Align
# Produce full outputs without version or citation information
# Save rotation matrices to ROT_DIR
# Save outputs to file
awk -v PDB_DIR="${PDB_DIR}" \
    -v ROT_DIR="${ROT_DIR}" \
    '{ system("./TMalign_cpp " PDB_DIR $1 " " PDB_DIR $2 " -outfmt -1 -m " ROT_DIR $1 "-" $2 ".txt") }' \
    $PAIR_FILE >> $OUT_FILE


echo "Completed: " $(date)