def param_log(){
if (params.gen_org=='human')
  log.info """
______________________________________________________

                ATAC PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
--concat_lanes                  ${params.concat_lanes}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--bowtie2Index                  ${params.bowtie2Index}
--bowtieMaxInsert               ${params.bowtieMaxInsert}
--bowtieVSensitive              ${params.bowtieVSensitive}
--cutadaptMinLength             ${params.cutadaptMinLength}
--cutadaptQualCutoff            ${params.cutadaptQualCutoff}
--cutadaptAdapterR1             ${params.cutadaptAdapterR1}
--cutadaptAdapterR2             ${params.cutadaptAdapterR2}
--tmpdir                        ${params.tmpdir}

Project Directory: ${projectDir}
______________________________________________________
"""
else
log.info """
______________________________________________________

                ATAC PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
--concat_lanes                  ${params.concat_lanes}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--effective_genome_size         ${params.effective_genome_size} 
--chain                         ${params.chain}
--bowtie2Index                  ${params.bowtie2Index}
--bowtieMaxInsert               ${params.bowtieMaxInsert}
--bowtieVSensitive              ${params.bowtieVSensitive}
--cutadaptMinLength             ${params.cutadaptMinLength}
--cutadaptQualCutoff            ${params.cutadaptQualCutoff}
--cutadaptAdapterR1             ${params.cutadaptAdapterR1}
--cutadaptAdapterR2             ${params.cutadaptAdapterR2}
--tmpdir                        ${params.tmpdir}

Project Directory: ${projectDir}
______________________________________________________
"""

}
