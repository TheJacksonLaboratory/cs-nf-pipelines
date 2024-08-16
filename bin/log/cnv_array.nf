def param_log() {
    // Check required parameters
    if (!params.bpm_file) {
        error "'--bpm_file': is not provided, it is a required parameter."
    }

    if (!params.egt_file) {
        error "'--egt_file': is not provided, it is a required parameter."
    }

    if (!params.gtc_csv) {
        error "'--gtc_csv': is not provided, it is a required parameter."
    }

    if (!params.ref_fa) {
        error "'--ref_fa': is not provided, it is a required parameter."
    }


    // Log parameter information
    log.info """
    CNV_ARRAY PARAMETER LOG

    --comment: ${params.comment ?: 'N/A'}

    Results Published to: ${params.pubdir ?: 'N/A'}
    ______________________________________________________
    --bpm_file                  ${params.bpm_file}
    --egt_file                  ${params.egt_file}
    --gtc_csv                   ${params.gtc_csv}
    --ref_fa                    ${params.ref_fa}
    -w                          ${workDir}
    --keep_intermediate         ${params.keep_intermediate ?: 'N/A'}
    -c                          ${params.config ?: 'N/A'}
    
    Project Directory: ${projectDir}
    Command line call: 
    ${workflow.commandLine}
    ______________________________________________________
    """
}
