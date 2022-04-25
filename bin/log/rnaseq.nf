def param_log(){
if (params.gen_org=='human')
  log.info """
______________________________________________________

                RNASEQ PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow             ${params.workflow}
--gen_org              ${params.gen_org}
--read_type            ${params.read_type}
--sample_folder        ${params.sample_folder}
--extension            ${params.extension}
--pattern              ${params.pattern}
--organize_by          ${params.organize_by}
--pubdir               ${params.pubdir}
-w                     ${workDir}
--keep_intermediate    ${params.keep_intermediate}
-c                     ${params.config}
--read_prep            ${params.read_prep}
--ref_fa               ${params.ref_fa}
--ref_fai              ${params.ref_fai}
--min_pct_hq_reads     ${params.min_pct_hq_reads}
--seed_length          ${params.seed_length}
--rsem_ref_prefix      ${params.rsem_ref_prefix}
--rsem_ref_files       ${params.rsem_ref_files}
--rsem_aligner         ${params.rsem_aligner}
--picard_dict          ${params.picard_dict}
--ref_flat             ${params.ref_flat}
--ribo_intervals       ${params.ribo_intervals}
--probes               ${params.probes}
--ctp_genes            ${params.ctp_genes}

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
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--read_prep    	                ${params.read_prep}
--ref_fa                        ${params.ref_fa}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--seed_length                   ${params.seed_length}
--rsem_ref_prefix               ${params.rsem_ref_prefix}
--rsem_ref_files                ${params.rsem_ref_files}
--rsem_aligner                  ${params.rsem_aligner}
--picard_dict                   ${params.picard_dict}

Project Directory: ${projectDir}
______________________________________________________
"""

}
