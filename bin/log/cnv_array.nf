import Logos

logo = new Logo()
println '\n'
println logo.show()

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
    --snp_platform              ${params.snp_platform}
    --gc_file                   ${params.gc_file}
    --rt_file                   ${params.rt_file}
    --chrArm                    ${params.chrArm}
    --cnvGeneFile               ${params.cnvGeneFile}
    -w                          ${workDir}
    -c                          ${params.config}

    Project Directory: ${projectDir}
    Command line call: 
    ${workflow.commandLine}
    ______________________________________________________
    """
}
