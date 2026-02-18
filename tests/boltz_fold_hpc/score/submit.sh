#!/bin/bash
#SBATCH --job-name="boltz_had_score"
##SBATCH --partition=Pisces
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=1:00:00


module load singularity
export SINGULARITY_CONTAINER_HOME=/scratch/cford38/haddock3_resources

singularity exec $SINGULARITY_CONTAINER_HOME/haddock3.sif python score.py --pdb_path $1
