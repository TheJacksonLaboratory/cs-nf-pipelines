import Logos

logo = new Logo()
println '\n'
println logo.show()

if (params.gen_org != "mouse" && params.gen_org != "human") {
  error "'--gen_org': \"${params.gen_org}\" is not valid, supported options are 'human' or 'mouse'" 
}

def param_log(){
if (params.gen_org=='human')
log.info """
PTA PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--csv_input                     ${params.csv_input}
--organize_by                   ${params.organize_by}
--pubdir                        ${params.pubdir}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--gen_org                       ${params.gen_org}
--pdx                           ${params.pdx}
--read_type                     ${params.read_type}
--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--detect_adapter_for_pe         ${params.detect_adapter_for_pe}
--xenome_prefix                 ${params.xenome_prefix}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--ref_fa_dict                   ${params.ref_fa_dict}
--combined_reference_set        ${params.combined_reference_set}
--mismatch_penalty              ${params.mismatch_penalty}
--gold_std_indels               ${params.gold_std_indels}
--phase1_1000G                  ${params.phase1_1000G}
--dbSNP                         ${params.dbSNP}
--dbSNP_index                   ${params.dbSNP_index}
--chrom_contigs                 ${params.chrom_contigs}
--chrom_intervals               ${params.chrom_intervals}
--call_val                      ${params.call_val}
--ploidy_val                    ${params.ploidy_val}
--excludeIntervalList           ${params.excludeIntervalList}
--hapmap                        ${params.hapmap}
--omni                          ${params.omni}
--pon_bed                       ${params.pon_bed}
--intervalListBed               ${params.intervalListBed}
--lancet_beds_directory         ${params.lancet_beds_directory}
--mappability_directory         ${params.mappability_directory}
--bicseq2_chromList             ${params.bicseq2_chromList}
--bicseq2_no_scaling            ${params.bicseq2_no_scaling}
--germline_filtering_vcf        ${params.germline_filtering_vcf}
--gripss_pon                    ${params.gripss_pon}
--callRegions                   ${params.callRegions}
--strelka_config                ${params.strelka_config}
--msisensor_model               ${params.msisensor_model}
--vep_cache_directory           ${params.vep_cache_directory}
--vep_fasta                     ${params.vep_fasta}
--cosmic_cgc                    ${params.cosmic_cgc}
--cosmic_cancer_resistance_muts ${params.cosmic_cancer_resistance_muts}
--ensembl_entrez                ${params.ensembl_entrez}
--cytoband                      ${params.cytoband}
--dgv                           ${params.dgv}
--thousandG                     ${params.thousandG}
--cosmicUniqueBed               ${params.cosmicUniqueBed}
--cancerCensusBed               ${params.cancerCensusBed}
--ensemblUniqueBed              ${params.ensemblUniqueBed}
--gap                           ${params.gap}
--dgvBedpe                      ${params.dgvBedpe}
--thousandGVcf                  ${params.thousandGVcf}
--svPon                         ${params.svPon}
--cosmicBedPe                   ${params.cosmicBedPe}
--na12878_bam                   ${params.na12878_bam}
--na12878_bai                   ${params.na12878_bai}
--na12878_sampleName            ${params.na12878_sampleName}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""
else
log.info """
PTA PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--csv_input                     ${params.csv_input}
--organize_by                   ${params.organize_by}
--pubdir                        ${params.pubdir}
-w                              ${workDir}
--keep_intermediate             ${params.keep_intermediate}
-c                              ${params.config}
--gen_org                       ${params.gen_org}
--read_type                     ${params.read_type}
--quality_phred                 ${params.quality_phred}
--unqualified_perc              ${params.unqualified_perc}
--detect_adapter_for_pe         ${params.detect_adapter_for_pe}
--xenome_prefix                 ${params.xenome_prefix}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--ref_fa_dict                   ${params.ref_fa_dict}
--combined_reference_set        ${params.combined_reference_set}
--mismatch_penalty              ${params.mismatch_penalty}
--dbSNP                         ${params.dbSNP}
--dbSNP_index                   ${params.dbSNP_index}
--chrom_contigs                 ${params.chrom_contigs}
--chrom_intervals               ${params.chrom_intervals}
--call_val                      ${params.call_val}
--ploidy_val                    ${params.ploidy_val}
--excludeIntervalList           ${params.excludeIntervalList}
--intervalListBed               ${params.intervalListBed}
--lancet_beds_directory         ${params.lancet_beds_directory}
--delly_exclusion               ${params.delly_exclusion}
--delly_mappability             ${params.delly_mappability}
--cnv_window                    ${params.cnv_window}
--cnv_min_size                  ${params.cnv_min_size}
--cnv_germline_prob             ${params.cnv_germline_prob}
--callRegions                   ${params.callRegions}
--strelka_config                ${params.strelka_config}
--vep_cache_directory           ${params.vep_cache_directory}
--vep_fasta                     ${params.vep_fasta}
--cytoband                      ${params.cytoband}
--known_del                     ${params.known_del}
--known_ins                     ${params.known_ins}
--known_inv                     ${params.known_inv}
--ensemblUniqueBed              ${params.ensemblUniqueBed}
--gap                           ${params.gap}
--exclude_list                  ${params.exclude_list}
--proxy_normal_bam              ${params.proxy_normal_bam}
--proxy_normal_bai              ${params.proxy_normal_bai}
--proxy_normal_sampleName       ${params.proxy_normal_sampleName}

Project Directory: ${projectDir}

Command line call: 
${workflow.commandLine}
______________________________________________________
"""

}
