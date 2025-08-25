
nextflow.preview.output = true

include { SEQ_LIST_TO_FASTA as MSA_SEQ_LIST_TO_FASTA; SEQ_LIST_TO_FASTA as INF_SEQ_LIST_TO_FASTA; NOOP_DEP } from './modules/tgen/af3'
include { MSA_WORKFLOW; INFERENCE_WORKFLOW } from './subworkflows/tgen/af3'
include { splitParquet } from 'plugin/nf-parquet'

def hash_from(String seq) {
    def md = java.security.MessageDigest.getInstance("SHA-256")
    md.update(seq.getBytes('UTF-8'))      // use getBytes(...) rather than .bytes
    def bytes = md.digest()
    bytes.collect { String.format('%02x', it) }.join()
}

process NOOP_DEP_META {

    label "af3_process_local"

    input:
    tuple val(meta), path(input_val)
    val value

    output:
    tuple val(meta), path("*", includeInputs: true)

    script:
    """
    true
    """
}


process CLEAN_BP3C50ID {
  label 'epident_local'

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
  label "epident_local"

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
  label 'bp3_inference'
  tag "bp3"
  
  input:
  path(fasta)
  
  output:
  path("*.csv")
  
  script:
  """
  export TORCH_HOME=${params.bp3_torch_home}
  
  bepipred3_CLI.py \\
      -i ${fasta} \\
      -o . \\
      -pred mjv_pred \\
      -add_seq_len \\
      -esm_dir ${params.bp3_esm_dir}
  """
}

process JOIN_BEPIPRED_INFERENCE {
    label "epident_local"

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

process EXTRACT_RSA {
    label "epident_local"

    input:
        path(parquet)
        val(inf_dir)
    
    output:
        path("*.parquet")
    
    script:
    """
    rsa_calculator.py \\
        --input_parquet ${parquet} \\
        --output_path ${parquet.getSimpleName()}.rsa.parquet \\
        --inference_path ${inf_dir}
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
    bp3_pq_bp3 = JOIN_BEPIPRED_INFERENCE.out
}

workflow AF3_INFERENCE_FROM_BP3_PARQUET {
    take:
    bp3_pq
    inf_dir

    main:

    af3_msa_ch = bp3_pq.splitParquet()
    .map{
        row -> 

        def seq_hash = hash_from(row["seq"])

        tuple(
            [
                id : seq_hash,
                protein_type : "any",
            ],
            [row["seq"]],
        )
    }

    af3_msa_fasta_ch = MSA_SEQ_LIST_TO_FASTA(af3_msa_ch)

    MSA_WORKFLOW(af3_msa_fasta_ch)

    new_meta_msa = MSA_WORKFLOW.out.new_meta_msa

    af3_inf_ch = bp3_pq.splitParquet()
    .map{
        row -> 

        tuple(
            [
                id : row["job_name"],
                protein_types : ["any"],
                segids : ["A"],
            ],
            [row["seq"]],
        )
    }

    af3_inf_fasta_ch = INF_SEQ_LIST_TO_FASTA(af3_inf_ch)

    af3_inf_fasta_ch = NOOP_DEP_META(af3_inf_fasta_ch, new_meta_msa.toList())

    INFERENCE_WORKFLOW(af3_inf_fasta_ch, inf_dir)

    emit:
    new_meta_inf = INFERENCE_WORKFLOW.out.new_meta_inf
}

workflow {
  
  main:
  test_fasta = Channel.fromPath(params.test_fasta)
  train_fasta = Channel.fromPath(params.train_fasta)

  inf_dir = file("$workflow.outputDir/inference").toUriString()

  CLEAN_BP3C50ID(train_fasta, test_fasta)

  bp3_pq = CLEAN_BP3C50ID.out.filt_pq
  discard_pq = CLEAN_BP3C50ID.out.discard_pq

  AF3_INFERENCE_FROM_BP3_PARQUET(bp3_pq, inf_dir)
  meta_inf = AF3_INFERENCE_FROM_BP3_PARQUET.out.new_meta_inf

  bp3_pq_af3 = NOOP_DEP(bp3_pq, meta_inf.toList())

  bp3_pq_rsa = EXTRACT_RSA(bp3_pq_af3, inf_dir)

  BP3_INFERENCE_FROM_BP3_PARQUET(bp3_pq)

  bp3_pq_bp3 = BP3_INFERENCE_FROM_BP3_PARQUET.out.bp3_pq_bp3

  publish:
    bp3_pq = bp3_pq
    discard_pq = discard_pq
    meta_inf = meta_inf
    bp3_pq_bp3 = bp3_pq_bp3
    bp3_pq_rsa = bp3_pq_rsa
}

output {
    bp3_pq {}
    discard_pq {}
    meta_inf { path "inference"}
    bp3_pq_bp3 {}
    bp3_pq_rsa {}
}