#!/bin/bash
#SBATCH --job-name=BP3C50ID
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/bp3c50id/slurm.%j.log

export NXF_LOG_FILE=logs/bp3c50id/.nextflow.log
export NXF_CACHE_DIR=logs/bp3c50id/.nextflow

conda run -n nf-core --live-stream nextflow run \
    ./workflows/bp3c50id.nf \
    --train_fasta data/bp3c50id/raw/BP3C50ID_training_set.fasta \
    --test_fasta data/bp3c50id/raw/BP3C50ID_external_test_set.fasta \
    -output-dir data/bp3c50id \
    -profile gemini \
    -resume