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
--non_directional      ${params.non_directional}
--trimLength           ${params.trimLength}
--qualThreshold        ${params.qualThreshold}
--adapOverlap          ${params.adapOverlap}
--adaptorSeq           ${params.adaptorSeq}
--seedLength           ${params.seedLength}
--seedMismatch         ${params.seedMismatch}
--MinInsert            ${params.MinInsert}
--MaxInsert            ${params.MaxInsert}
--ref_fa_index         ${params.ref_fa_index}
--aligner              ${params.aligner}
--skip_deduplication   ${params.skip_deduplication}
--cytosine_report      ${params.cytosine_report}
--comprehensive        ${params.comprehensive}

Project Directory: ${projectDir}
______________________________________________________
"""
else
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
--pubdir               ${params.pubdir}
--organize_by          ${params.organize_by}
-w                     ${workDir}
--keep_intermediate    ${params.keep_intermediate}
-c                     ${params.config}
--non_directional      ${params.non_directional}
--trimLength           ${params.trimLength}
--qualThreshold        ${params.qualThreshold}
--adapOverlap          ${params.adapOverlap}
--adaptorSeq           ${params.adaptorSeq}
--seedLength           ${params.seedLength}
--seedMismatch         ${params.seedMismatch}
--MinInsert            ${params.MinInsert}
--MaxInsert            ${params.MaxInsert}
--ref_fa_index         ${params.ref_fa_index}
--aligner              ${params.aligner}
--skip_deduplication   ${params.skip_deduplication}
--cytosine_report      ${params.cytosine_report}
--comprehensive        ${params.comprehensive}

Project Directory: ${projectDir}
______________________________________________________
"""

}
