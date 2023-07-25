import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
log.info """
GENERATE DO GBRS INPUT FILES

Results Published to: ${params.pubdir}
______________________________________________________
-w                          ${workDir}
--num_generations           ${params.num_generations}
--ensembl_build             ${params.ensembl_build}
--emase_gene2transcript     ${params.emase_gene2transcript}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
}