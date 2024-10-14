#!/bin/bash
#
#SBATCH --job-name=learnAlphabet
#
#SBATCH --time=8:00:00
#
#SBATCH -G 1
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --partition=jamesz
#
#SBATCH --output=slurm_outputs/%x.%j.out
#SBATCH --error=slurm_errors/%x.%j.out


echo "SLURM_JOBID" $SLURM_JOBID

source /home/users/shiye/.bashrc
conda activate strucsearch

cd /home/groups/jamesz/shiye/foldseek-analysis/training


# srun ./learnAlphabet.sh 20 100 data/pdbs_train_1000.txt data/pdbs_val.txt alphabet1/

# srun ./learnAlphabet.sh 20 25 data/pdbs_train.txt data/pdbs_val.txt alphabet3/

# srun ./learnAlphabet.sh 20 1 data/pdbs_train.txt data/pdbs_val.txt alphabet4/

srun ./learnAlphabet.sh 20 1 data/pdbs_train.txt data/pdbs_val.txt alphabet5/
