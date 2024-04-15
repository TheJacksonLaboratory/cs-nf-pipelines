import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
log.info """
WES PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--read_type                     ${params.read_type}
--csv_input                     ${params.csv_input}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--pdx                           ${params.pdx}
--xenome_index                  ${params.xenome_prefix}
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
--sequenza_gc                   ${params.sequenza_gc}
--ensembl_database              ${params.ensembl_database}
--contam_ref                    ${params.contam_ref}
--dbSNP                         ${params.dbSNP}
--target_gatk                   ${params.target_gatk}
--target_picard                 ${params.target_picard}
--bait_picard                   ${params.bait_picard}
--snpEff_config                 ${params.snpEff_config}
--mismatch_penalty              ${params.mismatch_penalty}
--call_val                      ${params.call_val}
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