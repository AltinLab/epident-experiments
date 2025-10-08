#!/bin/bash
#SBATCH --job-name=iedb_bp3
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/iedb_bp3/slurm.%j.log

export NXF_LOG_FILE=logs/iedb_bp3/.nextflow.log
export NXF_CACHE_DIR=logs/iedb_bp3/.nextflow

conda run -n nf-core --live-stream nextflow run \
    ./workflows/iedb_bp3.nf \
    --input data/iedb_bp3/raw \
    --inf_dir data/iedb_bp3/inference \
    -output-dir data/iedb_bp3 \
    -profile gemini \
    -resume