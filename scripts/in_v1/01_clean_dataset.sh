#!/bin/bash
#SBATCH --job-name=in_v1_clean_dataset
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/in_v1/clean_dataset/slurm.%j.log


conda run -n epident-experiments --live-stream python \
    ./workflows/bin/clean_in_v1.py \
    --fasta_path data/in_v1/raw/INF1_targets.fasta \
    --epitope_data_dir data/in_v1/raw \
    --output_path data/in_v1/unclustered