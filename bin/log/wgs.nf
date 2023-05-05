import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
if (params.gen_org=='human')
log.info """
WGS PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
________________________________________________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--gen_ver                       ${params.gen_ver}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
--concat_lanes                  ${params.concat_lanes}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--hq_pct                        ${params.hq_pct}
--dbSNP                         ${params.dbSNP}
--snpEff_config                 ${params.snpEff_config}
--mismatch_penalty              ${params.mismatch_penalty}
--call_val                      ${params.call_val}
--ploidy_val                    ${params.ploidy_val}
--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbNSFP                        ${params.dbNSFP}
--cosmic                        ${params.cosmic}
--snpEff_config                 ${params.snpEff_config}


Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
________________________________________________________________________________________
"""
else
log.info """
WGS PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
________________________________________________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--gen_ver                       ${params.gen_ver}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
--concat_lanes                  ${params.concat_lanes}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--hq_pct                        ${params.hq_pct}
--dbSNP                         ${params.dbSNP}
--snpEff_config                 ${params.snpEff_config}
--mismatch_penalty              ${params.mismatch_penalty}
--call_val                      ${params.call_val}
--ploidy_val                    ${params.ploidy_val}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
________________________________________________________________________________________
"""

}
