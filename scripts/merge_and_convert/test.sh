#!/bin/bash
#SBATCH --job-name="hf_merge_test"
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --gres=gpu:H200
#SBATCH --time=1:00:00

python test_model.py
