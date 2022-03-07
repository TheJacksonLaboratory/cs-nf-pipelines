def param_log(){
if (params.gen_org=='human')
  log.info """
______________________________________________________

                WES PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--filter_trim                   ${params.filter_trim}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--read_group_pyfile             ${params.read_group_pyfile}
--dbSNP                         ${params.dbSNP}
--target_gatk                   ${params.target_gatk}
--target_picard                 ${params.target_picard}
--bait_picard                   ${params.bait_picard}
--snpEff_config                 ${params.snpEff_config}
--stats_agg                     ${params.stats_agg}
--mismatch_penalty              ${params.mismatch_penalty}
--call_val                      ${params.call_val}
--ploidy_val                    ${params.ploidy_val}
--gen_ver                       ${params.gen_ver}
--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbNSFP                        ${params.dbNSFP}
--cosmic                        ${params.cosmic}
--cosmic_annot                  ${params.cosmic_annot}
--snpEff_config                 ${params.snpEff_config}


Project Directory: ${projectDir}
______________________________________________________
"""
else
log.info """
______________________________________________________

                WES PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--filter_trim                   ${params.filter_trim}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--read_group_pyfile             ${params.read_group_pyfile}
--dbSNP                         ${params.dbSNP}
--target_gatk                   ${params.target_gatk}
--target_picard                 ${params.target_picard}
--bait_picard                   ${params.bait_picard}
--snpEff_config                 ${params.snpEff_config}
--stats_agg                     ${params.stats_agg}
--mismatch_penalty              ${params.mismatch_penalty}
--call_val                      ${params.call_val}
--ploidy_val                    ${params.ploidy_val}

Project Directory: ${projectDir}
______________________________________________________
"""

}
