#!/bin/bash
#SBATCH --job-name=REPLACEME
#SBATCH --mail-type=ALL
#SBATCH --mail-user=REPLACEME
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=8:00:00
#SBATCH --output=tmp/nextflow/REPLACEME

export NXF_LOG_FILE=tmp/nextflow/REPLACEME
export NXF_CACHE_DIR=tmp/nextflow/REPLACEME

conda run -n nf-core --live-stream nextflow run \
    REPLACE_ME