import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
    log.info """
    GBRS RUN PARAMETER LOG

    --comment: ${params.comment}

    Results Published to: ${params.pubdir}
    ______________________________________________________
    --workflow                      ${params.workflow}
    -w                              ${workDir}
    -c                              ${params.config}
    --csv_input                     ${params.csv_input}
    --gen_org                       ${params.gen_org}
    --genome_build                  ${params.genome_build}
    --ref_fa                        ${params.ref_fa}
    --chrom_contigs                 ${params.chrom_contigs}

    Project Directory: ${projectDir}

    Command line call: 
    ${workflow.commandLine}
    ______________________________________________________
    """


}
