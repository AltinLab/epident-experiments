process SEQ_LIST_TO_FASTA {
    label "af3_process_local"

    input:
      tuple val(meta), val(seq_list)

    output:
      tuple val(meta), path("${meta.id}.fasta")

    script:
    """
    filename="${meta.id}.fasta"
    : > "\$filename"

    i=1
    for seq in ${seq_list.join(' ')}; do
        echo ">\$i"   >> "\$filename"
        echo "\$seq"  >> "\$filename"
        i=\$(( i + 1 ))
    done
    """
}

process NOOP_DEP {

    label "af3_process_local"

    input:
    path input_val
    val msa_done_key

    output:
    path "*", includeInputs: true

    script:
    """
    true
    """
}


process FILT_FORMAT_MSA {
    label "af3_process_local"
    tag "${meta.protein_type}-${meta.id}"
    conda "${moduleDir}/environment.yaml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json"), optional: true


    script:
    """
   filt_format_msa.py \\
        --protein_type "${meta.protein_type}" \\
        --job_name "${meta.id}" \\
        --fasta "$fasta" \\
        --msa_cache_dir "${params.msa_cache_dir}"
    """
}

process RUN_MSA {
    label "alphafold3_msa"

    memory { "${ Math.min(512, 64 * Math.pow(2, task.attempt - 1)) }GB" }
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    storeDir "${params.msa_cache_dir}/${meta.protein_type}"

    tag "${meta.protein_type}-${meta.id}"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("${meta.id}.json")

    script:
    """
    /app/alphafold/run_alphafold.py \\
        --json_path=$json \\
        --model_dir=${params.af3_model_dir} \\
        --db_dir=${params.af3_db_dir} \\
        --output_dir=. \\
        --norun_inference

    mv ${meta.id}/${meta.id}_data.json ${meta.id}.json
    """
}


process COMPOSE_INFERENCE_JSON {
    label "af3_process_local"
    conda "${moduleDir}/environment.yaml"
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta)
    // val because we dont want this to be resolved to relative
    val(inf_dir)

    output:
    tuple val(meta), path("*.json"), optional: true


    script:
    def segids = (meta.containsKey('segids')) ? "--segids ${meta.segids.join(',')}" : ''
    def skip_msa_arg = meta.containsKey("skip_msa") ? "--skip_msa ${meta.skip_msa.join(',')}" : ''
    """
    compose_inference_JSON.py \\
        --job_name "${meta.id}" \\
        --fasta "$fasta" \\
        --protein_types "${meta.protein_types.join(',')}" \\
        --msa_cache_dir "${params.msa_cache_dir}" \\
        --inf_dir "$inf_dir" \\
        ${segids} \\
        ${skip_msa_arg} \\
        --seeds 1,2,3,4,5 \\
        --check_inf_exists
    """
 }

process INFERENCE {
    tag "inference"
    label "alphafold3_inference"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("*")

    script:
    """
    python /app/alphafold/run_alphafold.py \\
        --json_path=$json \\
        --model_dir=$params.af3_model_dir \\
        --db_dir=$params.af3_db_dir \\
        --output_dir=. \\
        --norun_data_pipeline \\
        --save_embeddings \\
        --num_diffusion_samples=1
    """
}

process CLEAN_INFERENCE_DIR {
    label "af3_process_local"
    tag "clean_inference"
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5
    conda "${moduleDir}/environment.yaml"

    input:
    tuple val(meta), path(inference_dir)

    output:
    tuple val(meta), path("*", includeInputs: true)

    script:
    """
    clean_inference_dir.py \\
        -i $inference_dir
    """
}
