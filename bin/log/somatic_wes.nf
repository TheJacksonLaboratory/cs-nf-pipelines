import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
if (params.gen_org != "human") {
    error "'--gen_org': \"${params.gen_org}\" is not valid, supported option is 'human'" 
}
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
--pdx                           ${params.pdx}
--xengsort_host_fasta           ${params.xengsort_host_fasta}
--xengsort_idx_path             ${params.xengsort_idx_path}
--xengsort_idx_name             ${params.xengsort_idx_name}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--detect_adapter_for_pe         ${params.detect_adapter_for_pe}
--genotype_targets              ${params.genotype_targets}
--snpID_list                    ${params.snpID_list}
--snp_annotations               ${params.snp_annotations}
--snpweights_panel              ${params.snpweights_panel}
--gnomad_ref                    ${params.gnomad_ref}
--pon_ref                       ${params.pon_ref}
--genotype_pon                  ${params.genotype_pon}
--genotype_germline             ${params.genotype_germline}
--contam_ref                    ${params.contam_ref}
--dbSNP                         ${params.dbSNP}
--target_gatk                   ${params.target_gatk}
--target_picard                 ${params.target_picard}
--bait_picard                   ${params.bait_picard}
--snpEff_config                 ${params.snpEff_config}
--mismatch_penalty              ${params.mismatch_penalty}
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

}
