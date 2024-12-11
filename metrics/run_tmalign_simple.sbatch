#!/bin/bash
#
#SBATCH --job-name=run_tmalign_simple
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=jamesz,normal
#SBATCH --output=/home/groups/jamesz/shiye/foldseek-analysis/slurm_outputs/%x.%j.%a.out
#SBATCH --error=/home/groups/jamesz/shiye/foldseek-analysis/slurm_errors/%x.%j.%a.err
#SBATCH --time=5:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shiye@stanford.edu

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NAME="$SLURM_JOB_NAME
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory="$SLURM_SUBMIT_DIR

source /home/users/shiye/.bashrc

cd /home/groups/jamesz/shiye/foldseek-analysis

# command="./metrics/run_tmalign_simple.sh pairfile/pairfile_within_fold_train_all.out"
# command="./metrics/run_tmalign_simple.sh pairfile/pairfile_random_pairs_train_outside_fold_split_00.out"
# command="./metrics/run_tmalign_simple.sh pairfile/pairfile_random_pairs_train_outside_fold_split_01.out"
# command="./metrics/run_tmalign_simple.sh pairfile/pairfile_random_pairs_train_outside_fold_split_02.out"
# command="./metrics/run_tmalign_simple.sh pairfile/pairfile_random_pairs_train_outside_fold_split_03.out"
# command="./metrics/run_tmalign_simple.sh pairfile/pairfile_random_pairs_train_outside_fold_split_04.out"
# command="./metrics/run_tmalign_simple.sh pairfile/pairfile_random_pairs_train_outside_fold_split_05.out"
command="./metrics/run_tmalign_simple.sh pairfile/pairfile_within_fold_val_all.out"

echo "$command"
eval "$command"

exit 0
