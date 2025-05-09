import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
if (params.gen_org != "human") {
  error "'--gen_org': \"${params.gen_org}\" is not valid, supported option is 'human'" 
}

if (params.workflow=='amplicon_generic')

log.info """
AMPLICON PARAMETER LOG

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
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--multiqc_config                ${params.multiqc_config}

--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--detect_adapter_for_pe         ${params.detect_adapter_for_pe}

--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--mismatch_penalty              ${params.mismatch_penalty}

--markduplicates                ${params.markduplicates}

--amplicon_primer_intervals     ${params.amplicon_primer_intervals}
--amplicon_target_intervals     ${params.amplicon_target_intervals}

--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbSNP                         ${params.dbSNP}
--dbSNP_index                   ${params.dbSNP_index}
--ploidy_val                    ${params.ploidy_val}
--target_gatk                   ${params.target_gatk}
--call_val                      ${params.call_val}
--bwa_min_score                 ${params.bwa_min_score}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
else
"""
AMPLICON PARAMETER LOG

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
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--multiqc_config                ${params.multiqc_config}
--cutadaptMinLength             ${params.cutadaptMinLength}
--cutadaptQualCutoff            ${params.cutadaptQualCutoff}
--cutadaptAdapterR1             ${params.cutadaptAdapterR1}
--cutadaptAdapterR2             ${params.cutadaptAdapterR2}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--mismatch_penalty              ${params.mismatch_penalty}
--masterfile                    ${params.masterfile}
--amplicon_primer_intervals     ${params.amplicon_primer_intervals}
--amplicon_target_intervals     ${params.amplicon_target_intervals}
--amplicon_rsid_targets         ${params.amplicon_rsid_targets}
--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbSNP                         ${params.dbSNP}
--dbSNP_index                   ${params.dbSNP_index}
--ploidy_val                    ${params.ploidy_val}
--target_gatk                   ${params.target_gatk}
--call_val                      ${params.call_val}
--bwa_min_score                 ${params.bwa_min_score}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
}
