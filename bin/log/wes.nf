def param_log(){
if (params.gen_org=='human')
  log.info """
______________________________________________________

                RNASEQ PARAMITER LOG
______________________________________________________
--workflow                      ${params.workflow}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
-w                              ${workDir}
--cwd                           ${params.cwd}
--gen_org                       ${params.gen_org}
--extension                     ${params.extension}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--filter_trim                   ${params.filter_trim}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--read_group_pyfile             ${params.read_group_pyfile}
--stats_agg                     ${params.stats_agg}
--mismatch_penalty              ${params.mismatch_penalty}
--seed_length                   ${params.seed_length}
--rsem_ref_prefix               ${params.rsem_ref_prefix}
--rsem_ref_files                ${params.rsem_ref_files}
--rsem_aligner                  ${params.rsem_aligner}
--picard_dict                   ${params.picard_dict}
--summary_mets_PE               ${params.summary_mets_PE}
--summary_mets_SE               ${params.summary_mets_SE}
--probes                        ${params.probes}
--ctp_genes                     ${params.ctp_genes}
--gatk_form                     ${params.gatk_form}
--cov_calc                      ${params.cov_calc}
______________________________________________________
"""
else
log.info """
______________________________________________________

                RNASEQ PARAMITER LOG
______________________________________________________
--workflow                      ${params.workflow}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
-w                              ${workDir}
--cwd                           ${params.cwd}
--gen_org                       ${params.gen_org}
--extension                     ${params.extension}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--ref_fa                        ${params.ref_fa}
--filter_trim                   ${params.filter_trim}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--read_group_pyfile             ${params.read_group_pyfile}
--seed_length                   ${params.seed_length}
--rsem_ref_prefix               ${params.rsem_ref_prefix}
--rsem_ref_files                ${params.rsem_ref_files}
--rsem_aligner                  ${params.rsem_aligner}
--picard_dict                   ${params.picard_dict}
______________________________________________________
"""

}

