import Logos

logo = new Logo()
println '\n'
println logo.show()

def param_log(){
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
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--hq_pct                        ${params.hq_pct}
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

}
