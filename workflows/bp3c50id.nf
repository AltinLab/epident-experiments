
nextflow.preview.output = true

process CLEAN_BP3C50ID {
  label 'process_local'
  conda "envs/env.yaml"


  input:
    path train_fasta
    path test_fasta

  output:
    path("*discard*.parquet"), emit: discard_pq
    path("*filt*.parquet"), emit: filt_pq

  script:
  """
  clean_bp3c50id.py \\
    --train_fasta ${train_fasta} \\
    --test_fasta ${test_fasta} \\
    --discard_path bp3c50id.discard.parquet \\
    --output_path bp3c50id.filt.parquet
  """
}

process PARQUET_TO_FASTA {
  label "process_local"
  conda "envs/env.yaml"
  
  input:
      path(parquet)
  
  output:
      path("*.fasta")
  
  script:
  """
  #!/usr/bin/env python

  import polars as pl
  
  df = pl.read_parquet("${parquet}")
  
  with open("from_pq.fasta", "w") as f:
      for row in df.iter_rows(named=True):
          f.write(f">{row['job_name']}\\n{row['seq']}\\n")
  """
}

process RUN_BEPIPRED {
  queue 'gpu-v100'
  cpus '8'
  clusterOptions '--nodes=1 --ntasks=1 --gres=gpu:1 --time=2-00:00:00'
  memory '64GB'
  executor "slurm"
  tag "bp3"
  conda 'envs/bp3.yaml'
  
  input:
  path(fasta)
  
  output:
  path("*.csv")
  
  script:
  """
  export TORCH_HOME=${params.torch_home}
  
  bepipred3_CLI.py \\
      -i ${fasta} \\
      -o . \\
      -pred mjv_pred \\
      -add_seq_len \\
      -esm_dir ${params.esm_dir}
  """
}

process JOIN_BEPIPRED_INFERENCE {
    label "process_local"
    conda "envs/env.yaml"

    input:
        path(filt_dset)
        path(bepipred_csv)

    output:
        path("*.parquet")
    
    script:
    """
    #!/usr/bin/env python

    import polars as pl
    
    filt_dset = pl.read_parquet("${filt_dset}")
    bepipred_out = pl.read_csv("${bepipred_csv}").with_columns(
        pl.col("BepiPred-3.0 linear epitope score").str.strip_chars().cast(
            pl.Float64)
        )

    bp = (
        bepipred_out.group_by("Accession", maintain_order=True).agg(
            [
                pl.col("BepiPred-3.0 score").alias("bp3_score"),
                pl.col("BepiPred-3.0 linear epitope score").alias("bp3_linear_score"),
                # pl.col("Residue").str.join().alias("validate_seq"),
            ]
        )
    ).rename({"Accession": "job_name"})

    filt_dset = filt_dset.join(bp, on="job_name")
    filt_dset.write_parquet("${filt_dset.getSimpleName()}.bp3.parquet")
    """
}


workflow BP3_INFERENCE_FROM_BP3_PARQUET {
    take:
    bp3_pq

    main:
    PARQUET_TO_FASTA(bp3_pq)
    RUN_BEPIPRED(PARQUET_TO_FASTA.out)
    JOIN_BEPIPRED_INFERENCE(bp3_pq, RUN_BEPIPRED.out)

    emit:
    JOIN_BEPIPRED_INFERENCE.out
}

workflow AF3_INFERENCE_FROM_BP3_PARQUET {
    take:
    bp3_pq

    main:

    emit:

}

workflow {
  
  test_fasta = Channel.fromPath(params.test_fasta)
  train_fasta = Channel.fromPath(params.train_fasta)

  CLEAN_BP3C50ID(train_fasta, test_fasta)

  bp3_pq = CLEAN_BP3C50ID.filt_pq
  discard_pq = CLEAN_BP3C50ID.discard_pq

  output:
    bp3_pq = bp3_pq
    discard_pq = discard_pq

}

publish {
    bp3_pq
    discard_pq
}