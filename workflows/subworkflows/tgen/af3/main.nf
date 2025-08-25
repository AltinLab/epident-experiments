
include { FILT_FORMAT_MSA;
            RUN_MSA;
            COMPOSE_INFERENCE_JSON;
            INFERENCE;
            CLEAN_INFERENCE_DIR} from '../../../modules/tgen/af3'





workflow MSA_WORKFLOW {
    /*
    Arguments:
    - meta_fasta: Channel of tuples containing metadata and fasta files. Metadata should contain:
        - id: Unique identifier for the sequence
        - protein_types: List of protein types for each chain in the fasta (one of ["peptide", "mhc", "tcr", "any"])
        - segids (optional): List of segment IDs to label chains in the output structure
        - skip_msa (optional): List of integers to skip MSA for specific indices in the input fasta file

    Returns:
    - new_msa_list: List of newly completed MSAs. Can be used as a single token representing MSA completion
    */

    take:
    meta_fasta

    main:
    FILT_FORMAT_MSA(meta_fasta)
    RUN_MSA(FILT_FORMAT_MSA.out)

    emit:
    new_meta_msa = RUN_MSA.out

}

workflow INFERENCE_WORKFLOW {
    /*
    Arguments:
    - meta_fasta: Channel of tuples containing metadata and fasta files. Metadata should contain:
        - id: Unique identifier for the sequence
        - protein_types: List of protein types (one of ["peptide", "mhc", "tcr", "any"])
        - segids (optional): List of segment IDs to label chains in the output structure
    - inf_dir: Directory to check for existing inference results

    Returns:
    - new_inf_list: List of completed inferences. Can be used as a single token representing inference completion
    */

    take:
    meta_fasta
    inf_dir

    main:

    json = COMPOSE_INFERENCE_JSON(meta_fasta, inf_dir)

    inference = INFERENCE(json)


    inference = CLEAN_INFERENCE_DIR(inference)



    emit:
    new_meta_inf = inference

}


