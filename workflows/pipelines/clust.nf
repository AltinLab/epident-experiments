nextflow.preview.output = true

include { PARQUET_TO_FASTA;
            FASTA_TO_PARQUET;
            CLUSTER_FASTA;
            ANNOTATE_REPRESENTATIVES;} from '../modules/clust'

workflow {
    main:

    protein_pq = Channel.fromPath(params.input)

    PARQUET_TO_FASTA(protein_pq)
    CLUSTER_FASTA(PARQUET_TO_FASTA.out)
    FASTA_TO_PARQUET(CLUSTER_FASTA.out)
    ANNOTATE_REPRESENTATIVES(protein_pq, FASTA_TO_PARQUET.out)

    publish:
    representatives = ANNOTATE_REPRESENTATIVES.out
}

output {

    representatives {
        path "clustered"
    }
}