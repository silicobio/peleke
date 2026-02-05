#!/bin/bash
#SBATCH --job-name="hf_merge"
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --gres=gpu:H200
#SBATCH --time=1:00:00

# python merge_model.py --model_name silicobio/peleke-phi-4
# python merge_model.py --model_name silicobio/peleke-llama-3.1-8b-instruct
python merge_model.py --model_name silicobio/peleke-mistral-7b-instruct-v0.2
