def param_log() {
    // Check required parameters
    if (!params.bpm_file) {
        error "'--bpm_file': is not provided, it is a required parameter."
    }

    if (!params.egt_file) {
        error "'--egt_file': is not provided, it is a required parameter."
    }

    if (!params.csv_file) {
        error "'--csv_file': is not provided, it is a required parameter."
    }

    if (!params.fasta_file) {
        error "'--fasta_file': is not provided, it is a required parameter."
    }

    if (!params.tsv_file) {
        error "'--tsv_file': is not provided, it is a required parameter."
    }

    // Log parameter information
    log.info """
    CNV_ARRAY PARAMETER LOG

    --comment: ${params.comment ?: 'N/A'}

    Results Published to: ${params.pubdir ?: 'N/A'}
    ______________________________________________________
    --idat_folder               ${params.idat_folder ?: 'N/A'}
    --bpm_file                  ${params.bpm_file}
    --egt_file                  ${params.egt_file}
    --csv_file                  ${params.csv_file}
    --fasta_file                ${params.fasta_file}
    --tsv_file                  ${params.tsv_file}
    -w                          ${workDir}
    --keep_intermediate         ${params.keep_intermediate ?: 'N/A'}
    -c                          ${params.config ?: 'N/A'}
    
    Project Directory: ${projectDir}
    Command line call: 
    ${workflow.commandLine}
    ______________________________________________________
    """
}
