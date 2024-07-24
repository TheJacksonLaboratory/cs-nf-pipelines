import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){

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

log.info """
IAAP_CLI PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--csv_input                 ${params.csv_input}
--bpm_file                  ${params.bpm_file}
--egt_file                  ${params.egt_file}
-w                          ${workDir}
--keep_intermediate         ${params.keep_intermediate}
-c                          ${params.config}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
}


log.info """
BCFTOOLS_GTC2VCF PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--bpm_file                  ${params.bpm_file}
--csv_file                  ${params.csv_file}
--egt_file                  ${params.egt_file}
--fasta_file                ${params.fasta_file}
--tsv_file                  ${params.tsv_file}
-w                          ${workDir}
--keep_intermediate         ${params.keep_intermediate}
-c                          ${params.config}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
}
