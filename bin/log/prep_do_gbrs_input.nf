def param_log(){
log.info """
______________________________________________________

        GENERATE DO GBRS INPUT FILES

Results Published to: ${params.pubdir}
______________________________________________________
-w                          ${workDir}
--num_generations           ${params.num_generations}
--ensembl_build             ${params.ensembl_build}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
}