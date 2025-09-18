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
  
  with open("${parquet.getSimpleName()}.fasta", "w") as f:
      for row in df.iter_rows(named=True):
          f.write(f">{row['job_name']}\\n{row['seq']}\\n")
  """
}

process FASTA_TO_PARQUET {
  label "epident_local"
  
  input:
      path(fasta)
  
  output:
      path("*.parquet")
  
  script:
  """
  #!/usr/bin/env python

  import polars as pl
  from epident.utils import fasta_to_polars
  
  df = fasta_to_polars("${fasta}").rename({"name" : "job_name"})
  
  df.write_parquet("${fasta.getSimpleName()}.mmseqs.parquet")
  """
}


process CLUSTER_FASTA {
    label "mmseqs"

    input:
        path(fasta)

    output:
        path("*.fasta")

    script:
    """
    mmseqs createdb ${fasta} DB && \\
    mmseqs cluster DB DB_clu tmp \\
        --min-seq-id 0.7 && \\
    mmseqs createsubdb DB_clu DB DB_clu_rep && \\
    mmseqs convert2fasta DB_clu_rep "${fasta.getSimpleName()}.clust.fasta"
    """
}


process ANNOTATE_REPRESENTATIVES {
    label "epident_local"

    input:
        path(orig_pq)
        path(rep_pq)

    output:
        path("*.parquet")

    script:
    """
    #!/usr/bin/env python

    import polars as pl

    orig_pq = pl.read_parquet("${orig_pq}")

    rep_pq = pl.read_parquet("${rep_pq}").with_columns(pl.lit(True).alias("representative"))

    out_df = orig_pq.join(rep_pq.select(pl.exclude("seq")), on="job_name", how="left").with_columns(
        pl.when(pl.col("representative").is_null())
        .then(pl.lit(False))
        .otherwise(pl.col("representative"))
        .alias("representative")
    ).filter(pl.col("representative")).select(pl.exclude("representative"))

    out_df.write_parquet("${orig_pq.getSimpleName()}.clust.parquet")
    """
}