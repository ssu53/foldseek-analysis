#!/bin/bash -x

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"


# Set input arguments and defaults

ALIGNMENTS_DIR="${1:-$'../training/tmp19_v4_encodings/alignments'}"
NUM_SAMPLES="${2:-1000}"


# Set paths

PDB_DIR='/scratch/groups/jamesz/shiye/scope40/'
MASTER_FILE='MASTER_RESULTS.csv'
ROT_DIR='rot_mats/'
RAW_OUT_FILE='tmalign.out'
PARSED_OUT_FILE='tmalign.csv'
AUGMENTED_OUT_FILE='tmalign_extra.csv'
PAIR_FILE='pairfile.out'


# Put output artefacts in directory named by the datetime

DIR="outputs_$(date '+%Y-%m-%d-%H-%M-%S')"
ROT_DIR="${DIR}/${ROT_DIR}"
RAW_OUT_FILE="${DIR}/${RAW_OUT_FILE}"
PARSED_OUT_FILE="${DIR}/${PARSED_OUT_FILE}"
AUGMENTED_OUT_FILE="${DIR}/${AUGMENTED_OUT_FILE}"
PAIR_FILE="${DIR}/${PAIR_FILE}"

mkdir $DIR
mkdir $ROT_DIR


# Get the pairfile

# Run with an existing pair file instead of sampling it from alignments
# EXISTING_PAIR_FILE='/home/groups/jamesz/shiye/foldseek-analysis/training/data/tmaln-06_500.out'
# EXISTING_PAIR_FILE='/home/groups/jamesz/shiye/foldseek-analysis/metrics/pairfile_random_subset.out'
# EXISTING_PAIR_FILE='/home/groups/jamesz/shiye/foldseek-analysis/metrics/pairfile_train.out'
# EXISTING_PAIR_FILE='/home/groups/jamesz/shiye/foldseek-analysis/metrics/pairfile_val.out'
EXISTING_PAIR_FILE='/home/groups/jamesz/shiye/foldseek-analysis/metrics/pairfile_random_pairs_val.out'


if [ -z ${EXISTING_PAIR_FILE+x} ]; then 
    echo "Sampling pair file from " $ALIGNMENTS_DIR
    python get_sampled_pairfile.py $ALIGNMENTS_DIR $NUM_SAMPLES $MASTER_FILE $PAIR_FILE
else 
    echo "Using existing pair file " $EXISTING_PAIR_FILE
    cp $EXISTING_PAIR_FILE $PAIR_FILE
fi



# Log stuff to output file

echo "### Running TM-align (Version 20220412), compiled locally from cpp." > $RAW_OUT_FILE

echo "### PDB_DIR" $PDB_DIR >> $RAW_OUT_FILE
echo "### MASTER_FILE" $MASTER_FILE >> $RAW_OUT_FILE
echo "### ALIGNMENTS_DIR" $ALIGNMENTS_DIR >> $RAW_OUT_FILE
echo "### PAIR_FILE" $PAIR_FILE >> $RAW_OUT_FILE
echo "### ROT_DIR" $ROT_DIR >> $RAW_OUT_FILE
echo "### RAW_OUT_FILE" $RAW_OUT_FILE >> $RAW_OUT_FILE
echo "### PARSED_OUT_FILE" $PARSED_OUT_FILE >> $RAW_OUT_FILE
echo "### AUGMENTED_OUT_FILE" $AUGMENTED_OUT_FILE >> $RAW_OUT_FILE
echo "### PAIR_FILE" $PAIR_FILE >> $RAW_OUT_FILE

echo "Running..." $(date)



# Run TM-Align

# Produce full outputs without version or citation information
# Save rotation matrices to ROT_DIR
# Save outputs to file

# awk -v PDB_DIR="${PDB_DIR}" -v ROT_DIR="${ROT_DIR}" '{ system("echo " PDB_DIR $1 " " PDB_DIR $2 " -m " ROT_DIR $1 "-" $2 ".txt") }' $PAIR_FILE

awk -v PDB_DIR="${PDB_DIR}" \
    -v ROT_DIR="${ROT_DIR}" \
    '{ system("./TMalign_cpp " PDB_DIR $1 " " PDB_DIR $2 " -outfmt -1 -m " ROT_DIR $1 "-" $2 ".txt") }' \
    $PAIR_FILE >> $RAW_OUT_FILE

echo "Finished running TMalign: " $(date)



# Parse the TMalign outputs

python parse_tmalign_output.py $RAW_OUT_FILE $PARSED_OUT_FILE
echo "Parsed TMalign outputs: " $(date)



# Compute the LDDT, Chamfer, EMD metrics

python compute_extra_metrics.py $PARSED_OUT_FILE $PDB_DIR $ROT_DIR $AUGMENTED_OUT_FILE
echo "Computed extra metrics: " $(date)


# Update master results file by concatenating new results

awk FNR!=1 $AUGMENTED_OUT_FILE >> $MASTER_FILE
echo "Updated master results file: " $(date)
