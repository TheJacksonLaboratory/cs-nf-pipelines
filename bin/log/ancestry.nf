import Logos

def param_log(){

logo = new Logo()
println '\n'
println logo.show()

log.info """
WES PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--sample_folder                 ${params.sample_folder}
--csv_input                     ${params.csv_input}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--ref_fa                        ${params.ref_fa}
--genotype_targets              ${params.genotype_targets}
--snpID_list                    ${params.snpID_list}
--snp_annotations               ${params.snp_annotations}
--snpweights_panel              ${params.snpweights_panel}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

}
