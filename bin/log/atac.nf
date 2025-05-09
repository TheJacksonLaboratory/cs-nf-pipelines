import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
if (params.gen_org != "human" && params.gen_org != "mouse") {
  error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'mouse', or 'human'" 
}

if (params.gen_org=='human')
log.info """
ATAC PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
--concat_lanes                  ${params.concat_lanes}
--csv_input                     ${params.csv_input}
--download_data                 ${params.download_data}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--merge_replicates              ${params.merge_replicates}
--bowtie2Index                  ${params.bowtie2Index}
--bowtieMaxInsert               ${params.bowtieMaxInsert}
--bowtieVSensitive              ${params.bowtieVSensitive}
--cutadaptMinLength             ${params.cutadaptMinLength}
--cutadaptQualCutoff            ${params.cutadaptQualCutoff}
--cutadaptAdapterR1             ${params.cutadaptAdapterR1}
--cutadaptAdapterR2             ${params.cutadaptAdapterR2}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
else
log.info """
ATAC PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--genome_build                  ${params.genome_build}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
--concat_lanes                  ${params.concat_lanes}
--csv_input                     ${params.csv_input}
--download_data                 ${params.download_data}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--merge_replicates              ${params.merge_replicates}
--effective_genome_size         ${params.effective_genome_size} 
--chain                         ${params.chain}
--bowtie2Index                  ${params.bowtie2Index}
--bowtieMaxInsert               ${params.bowtieMaxInsert}
--bowtieVSensitive              ${params.bowtieVSensitive}
--cutadaptMinLength             ${params.cutadaptMinLength}
--cutadaptQualCutoff            ${params.cutadaptQualCutoff}
--cutadaptAdapterR1             ${params.cutadaptAdapterR1}
--cutadaptAdapterR2             ${params.cutadaptAdapterR2}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

}
