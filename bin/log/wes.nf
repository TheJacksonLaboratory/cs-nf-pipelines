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
WES PARAMETER LOG

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
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--hq_pct                        ${params.hq_pct}
--run_gvcf                      ${params.run_gvcf}
--dbSNP                         ${params.dbSNP}
--target_gatk                   ${params.target_gatk}
--target_picard                 ${params.target_picard}
--bait_picard                   ${params.bait_picard}
--snpEff_config                 ${params.snpEff_config}
--mismatch_penalty              ${params.mismatch_penalty}
--call_val                      ${params.call_val}
--ploidy_val                    ${params.ploidy_val}
--gen_ver                       ${params.gen_ver}
--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbNSFP                        ${params.dbNSFP}
--cosmic                        ${params.cosmic}
--snpEff_config                 ${params.snpEff_config}


Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
else
log.info """
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
--concat_lanes                  ${params.concat_lanes}
--csv_input                     ${params.csv_input}
--download_data                 ${params.download_data}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--hq_pct                        ${params.hq_pct}
--run_gvcf                      ${params.run_gvcf}
--dbSNP                         ${params.dbSNP}
--target_gatk                   ${params.target_gatk}
--target_picard                 ${params.target_picard}
--bait_picard                   ${params.bait_picard}
--snpEff_config                 ${params.snpEff_config}
--mismatch_penalty              ${params.mismatch_penalty}
--call_val                      ${params.call_val}
--ploidy_val                    ${params.ploidy_val}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

}
