#!/bin/bash
#SBATCH --job-name=hv_v1_clean_dataset
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/hv_v1/clean_dataset/slurm.%j.log


conda run -n epident-experiments --live-stream python \
    ./workflows/bin/clean_hv_v1.py \
    --fasta_path data/hv_v1/raw/fulldesign_2019-02-27_wGBKsw.fasta \
    --epitope_data_dir data/hv_v1/raw \
    --output_path data/hv_v1/unclustered