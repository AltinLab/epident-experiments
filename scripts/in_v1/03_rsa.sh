#!/bin/bash
#SBATCH --job-name=in_v1_rsa
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/in_v1/rsa/slurm.%j.log


conda run -n epident-experiments --live-stream python \
    ./workflows/bin/rsa_calculator.py \
    --input_parquet data/in_v1/clustered/in_v1_proteins.clust.parquet \
    --inference_path data/in_v1/inference \
    --output_path data/in_v1/clustered/in_v1_proteins.rsa.parquet