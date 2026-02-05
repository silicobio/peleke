#!/bin/bash
#SBATCH --job-name="boltz_peleke"
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --gres=gpu:1

##SBATCH --nodes=1
##SBATCH --ntasks-per-node=32
##SBATCH --mem=256GB
##SBATCH --gres=gpu:L40S:2
#SBATCH --time=1:00:00

module load singularity
# export SINGULARITY_CONTAINER_HOME=$(pwd)
export SINGULARITY_CONTAINER_HOME=/scratch/cford38/boltz2_resources
export BOLTZ_CACHE=/scratch/cford38/boltz2_resources

module load cuda/12.8

singularity exec --nv $SINGULARITY_CONTAINER_HOME/boltz2.sif boltz predict $1 --use_msa_server --accelerator gpu