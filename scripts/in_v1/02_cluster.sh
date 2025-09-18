#!/bin/bash
#SBATCH --job-name=in_v1_clust
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/in_v1/clust/slurm.%j.log

export NXF_LOG_FILE=logs/in_v1/clust/.nextflow.log
export NXF_CACHE_DIR=logs/in_v1/clust/.nextflow

conda run -n nf-core --live-stream nextflow run \
    ./workflows/pipelines/clust.nf \
    --input data/in_v1/unclustered/in_v1_proteins.parquet \
    -output-dir data/in_v1 \
    -profile gemini