def param_log(){
if (params.gen_org=='human')
  log.info """
______________________________________________________

                RRBS PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--gen_org              ${params.gen_org}
--read_type            ${params.read_type}
--sample_folder        ${params.sample_folder}
--extension            ${params.extension}
--pattern              ${params.pattern}
--concat_lanes         ${params.concat_lanes}
--organize_by          ${params.organize_by}
--pubdir               ${params.pubdir}
-w                     ${workDir}
--keep_intermediate    ${params.keep_intermediate}
-c                     ${params.config}


Project Directory: ${projectDir}
______________________________________________________
"""
else
log.info """
______________________________________________________

                RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--extension                     ${params.extension}
--pattern                       ${params.pattern}
--concat_lanes                  ${params.concat_lanes}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}


Project Directory: ${projectDir}
______________________________________________________
"""

}
