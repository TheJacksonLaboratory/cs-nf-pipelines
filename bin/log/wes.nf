def param_log(){
if (params.gen_org=='human')
  log.info """
______________________________________________________

                WES PARAMETER LOG
______________________________________________________
--workflow                      ${params.workflow}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--workDir                       ${workDir}
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

______________________________________________________
"""
else
log.info """
______________________________________________________

                WES PARAMETER LOG
______________________________________________________
--workflow                      ${params.workflow}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--w                             ${workDir}
--cwd                           ${params.cwd}
--gen_org                       ${params.gen_org}
--extension                     ${params.extension}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--ref_fa                        ${params.ref_fa}
--filter_trim                   ${params.filter_trim}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--read_group_pyfile             ${params.read_group_pyfile}

______________________________________________________
"""

}

