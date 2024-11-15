#!/bin/bash
#
#SBATCH --job-name=run_tmalign
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --partition=jamesz
#SBATCH --output=/home/groups/jamesz/shiye/foldseek-analysis/slurm_outputs/%x.%j.%a.out
#SBATCH --error=/home/groups/jamesz/shiye/foldseek-analysis/slurm_errors/%x.%j.%a.err
#SBATCH --time=3:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shiye@stanford.edu

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NAME="$SLURM_JOB_NAME
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory="$SLURM_SUBMIT_DIR

source /home/users/shiye/.bashrc

conda activate strucsearch

cd /home/groups/jamesz/shiye/foldseek-analysis

command="./metrics/run_tmalign.sh"

echo "$command"
eval "$command"

exit 0
