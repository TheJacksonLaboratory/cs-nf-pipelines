def param_log(){
if (params.gen_org=='human')
  log.info """
______________________________________________________

                RNASEQ PARAMETER LOG

--comment: ${params.comment}

Rsults Published to: ${params.pubdir}
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
-c                     ${params.config}
--ref_fa               ${params.ref_fa}
--ref_fai              ${params.ref_fai}
--ref_fa_indices       ${params.ref_fa_indices}
--filter_trim          ${params.filter_trim}
--min_pct_hq_reads     ${params.min_pct_hq_reads}
--read_group_pyfile    ${params.read_group_pyfile}
--stats_agg            ${params.stats_agg}
--mismatch_penalty     ${params.mismatch_penalty}
--seed_length          ${params.seed_length}
--rsem_ref_prefix      ${params.rsem_ref_prefix}
--rsem_ref_files       ${params.rsem_ref_files}
--rsem_aligner         ${params.rsem_aligner}
--picard_dict          ${params.picard_dict}
--ref_flat             ${params.ref_flat}
--ribo_intervals       ${params.ribo_intervals}
--summary_mets_PE      ${params.summary_mets_PE}
--summary_mets_SE      ${params.summary_mets_SE}
--probes               ${params.probes}
--ctp_genes            ${params.ctp_genes}
--gatk_form            ${params.gatk_form}
--cov_calc             ${params.cov_calc}

Project Directory: ${projectDir}
______________________________________________________
"""
else
log.info """
______________________________________________________

                RNASEQ PARAMETER LOG

--comment: ${params.comment}

Rsults Published to: ${params.pubdir}
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
-c                              ${params.config}
--ref_fa                        ${params.ref_fa}
--filter_trim                   ${params.filter_trim}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--read_group_pyfile             ${params.read_group_pyfile}
--seed_length                   ${params.seed_length}
--rsem_ref_prefix               ${params.rsem_ref_prefix}
--rsem_ref_files                ${params.rsem_ref_files}
--rsem_aligner                  ${params.rsem_aligner}
--picard_dict                   ${params.picard_dict}

Project Directory: ${projectDir}
______________________________________________________
"""

}
