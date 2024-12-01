#!/bin/bash -x


# Go to the correct directory

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"


# Pairfile to run

EXISTING_PAIR_FILE="${1:-$'pairfile/pairfile_within_fold_train_all.out'}"


# Set paths

PDB_DIR='/scratch/groups/jamesz/shiye/scope40/'
RAW_OUT_FILE='tmalign.out'
PARSED_OUT_FILE='tmalign.csv'
PAIR_FILE='pairfile.out'


# Put output artefacts in directory named by the datetime

DIR="outputs_$(date '+%Y-%m-%d-%H-%M-%S-%5N')"
RAW_OUT_FILE="${DIR}/${RAW_OUT_FILE}"
PARSED_OUT_FILE="${DIR}/${PARSED_OUT_FILE}"
PAIR_FILE="${DIR}/${PAIR_FILE}"

mkdir $DIR
cp $EXISTING_PAIR_FILE $PAIR_FILE


# Log stuff to output file

echo "### Running TM-align (Version 20220412), compiled locally from cpp." > $RAW_OUT_FILE

echo "### PDB_DIR" $PDB_DIR >> $RAW_OUT_FILE
echo "### RAW_OUT_FILE" $RAW_OUT_FILE >> $RAW_OUT_FILE
echo "### PARSED_OUT_FILE" $PARSED_OUT_FILE >> $RAW_OUT_FILE
echo "### PAIR_FILE" $PAIR_FILE >> $RAW_OUT_FILE


echo "Running..." $(date)

# Run TM-Align

# Produce compact tabular format output

awk -v PDB_DIR="${PDB_DIR}" \
    '{ system("./TMalign_cpp " PDB_DIR $1 " " PDB_DIR $2 " -outfmt 2") }' \
    $PAIR_FILE >> $RAW_OUT_FILE

echo "Finished running TMalign: " $(date)
