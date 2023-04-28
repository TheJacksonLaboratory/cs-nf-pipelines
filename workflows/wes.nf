#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/wes.nf"
include {param_log} from "${projectDir}/bin/log/wes.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {JAX_TRIMMER} from "${projectDir}/modules/utility_modules/jax_trimmer"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {GATK_BASERECALIBRATOR} from "${projectDir}/modules/gatk/gatk_baserecalibrator"
include {GATK_APPLYBQSR} from "${projectDir}/modules/gatk/gatk_applybqsr"
include {PICARD_COLLECTHSMETRICS} from "${projectDir}/modules/picard/picard_collecthsmetrics"
include {GATK_HAPLOTYPECALLER;
         GATK_HAPLOTYPECALLER as GATK_HAPLOTYPECALLER_GVCF} from "${projectDir}/modules/gatk/gatk_haplotypecaller"
include {GATK_VARIANTFILTRATION;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from "${projectDir}/modules/gatk/gatk_variantfiltration"
include {GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from "${projectDir}/modules/gatk/gatk_selectvariants"
include {GATK_MERGEVCF as GATK_MERGEVCF_UNANNOTATED;
         GATK_MERGEVCF as GATK_MERGEVCF_ANNOTATED} from "${projectDir}/modules/gatk/gatk_mergevcf"
include {GATK_INDEXFEATUREFILE} from "${projectDir}/modules/gatk/gatk_indexfeaturefile"
include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_DBSNP;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_DBSNP;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_DBSNP} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {SNPEFF;
         SNPEFF as SNPEFF_SNP;
         SNPEFF as SNPEFF_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_snpeff"
include {SNPEFF_ONEPERLINE;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_SNP;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_oneperline"
include {SNPSIFT_EXTRACTFIELDS} from "${projectDir}/modules/snpeff_snpsift/snpsift_extractfields"
include {SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpsift_dbnsfp"
include {AGGREGATE_STATS} from "${projectDir}/modules/utility_modules/aggregate_stats_wes"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"


// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// prepare reads channel
if (params.concat_lanes){
  if (params.read_type == 'PE'){
    read_ch = Channel
            .fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true, flat:true )
            .map { file, file1, file2 -> tuple(getLibraryId(file), file1, file2) }
            .groupTuple()
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}", checkExists:true, size:1 )
                .map { file, file1 -> tuple(getLibraryId(file), file1) }
                .groupTuple()
                .map{t-> [t[0], t[1].flatten()]}
  }
} else {
  if (params.read_type == 'PE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
  }
}

// if channel is empty give error message and exit
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern} and file extension: ${params.extension}"}

// main workflow
workflow WES {
  // Step 0: Concatenate Fastq files if required. 
  if (params.concat_lanes){
    if (params.read_type == 'PE'){
        CONCATENATE_READS_PE(read_ch)
        read_ch = CONCATENATE_READS_PE.out.concat_fastq
    } else if (params.read_type == 'SE'){
        CONCATENATE_READS_SE(read_ch)
        read_ch = CONCATENATE_READS_SE.out.concat_fastq
    }
  }

  // Step 1: Qual_Stat
  JAX_TRIMMER(read_ch)

  FASTQC(JAX_TRIMMER.out.trimmed_fastq)

  // Step 2: Get Read Group Information
  READ_GROUPS(JAX_TRIMMER.out.trimmed_fastq, "gatk")

  // Step 3: BWA-MEM Alignment
  if (params.read_type == 'PE') {
    bwa_mem_mapping = JAX_TRIMMER.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
  } else {
    bwa_mem_mapping = JAX_TRIMMER.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
  }

  BWA_MEM(bwa_mem_mapping)

  // Step 4: Variant Preprocessing - Part 1
  PICARD_SORTSAM(BWA_MEM.out.sam)
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)

  // If Human: Step 5-10
  ch_GATK_BASERECALIBRATOR_multiqc = Channel.empty() //optional log file for human only.
  if (params.gen_org=='human'){

    // Step 5: Variant Pre-Processing - Part 2
      GATK_BASERECALIBRATOR(PICARD_MARKDUPLICATES.out.dedup_bam)
      ch_GATK_BASERECALIBRATOR_multiqc = GATK_BASERECALIBRATOR.out.table // set log file for multiqc

      apply_bqsr = PICARD_MARKDUPLICATES.out.dedup_bam.join(GATK_BASERECALIBRATOR.out.table)
      GATK_APPLYBQSR(apply_bqsr)

    // Step 6: Variant Pre-Processing - Part 3
      collect_metrics = GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai)
      PICARD_COLLECTHSMETRICS(collect_metrics)

    // Step 7: Variant Calling
      haplotype_caller = GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai)
      GATK_HAPLOTYPECALLER(haplotype_caller, 'variant')

      haplotype_caller_gvcf = GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai)
      GATK_HAPLOTYPECALLER_GVCF(haplotype_caller_gvcf, 'gvcf')

    // Step 8: Variant Filtration
      // SNP
        select_var_snp = GATK_HAPLOTYPECALLER.out.vcf.join(GATK_HAPLOTYPECALLER.out.idx)
        GATK_SELECTVARIANTS_SNP(select_var_snp, 'SNP', 'selected_SNP')

        var_filter_snp = GATK_SELECTVARIANTS_SNP.out.vcf.join(GATK_SELECTVARIANTS_SNP.out.idx)
        GATK_VARIANTFILTRATION_SNP(var_filter_snp, 'SNP')

      // INDEL
        select_var_indel = GATK_HAPLOTYPECALLER.out.vcf.join(GATK_HAPLOTYPECALLER.out.idx)
      	GATK_SELECTVARIANTS_INDEL(select_var_indel, 'INDEL', 'selected_INDEL')

        var_filter_indel = GATK_SELECTVARIANTS_INDEL.out.vcf.join(GATK_SELECTVARIANTS_INDEL.out.idx)
        GATK_VARIANTFILTRATION_INDEL(var_filter_indel, 'INDEL')

    // Step 9: Post Variant Calling Processing - Part 1
      // SNP
        SNPSIFT_ANNOTATE_SNP_DBSNP(GATK_VARIANTFILTRATION_SNP.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')
        SNPSIFT_ANNOTATE_SNP_COSMIC(SNPSIFT_ANNOTATE_SNP_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmicID')
        SNPEFF_SNP(SNPSIFT_ANNOTATE_SNP_COSMIC.out.vcf, 'SNP', 'vcf')
        SNPSIFT_DBNSFP_SNP(SNPEFF_SNP.out.vcf, 'SNP')
        SNPEFF_ONEPERLINE_SNP(SNPSIFT_DBNSFP_SNP.out.vcf, 'SNP')

      // INDEL
        SNPSIFT_ANNOTATE_INDEL_DBSNP(GATK_VARIANTFILTRATION_INDEL.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')
        SNPSIFT_ANNOTATE_INDEL_COSMIC(SNPSIFT_ANNOTATE_INDEL_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmicID')
        SNPEFF_INDEL(SNPSIFT_ANNOTATE_INDEL_COSMIC.out.vcf, 'INDEL', 'vcf')
        SNPSIFT_DBNSFP_INDEL(SNPEFF_INDEL.out.vcf, 'INDEL')
        SNPEFF_ONEPERLINE_INDEL(SNPSIFT_DBNSFP_INDEL.out.vcf, 'INDEL')

    // Step 10: Post Variant Calling Processing - Part 2
      vcf_files_unannotated = SNPSIFT_ANNOTATE_SNP_COSMIC.out.vcf.join(SNPSIFT_ANNOTATE_INDEL_COSMIC.out.vcf)
      GATK_MERGEVCF_UNANNOTATED (vcf_files_unannotated, 'SNP_INDEL_filtered_unannotated_final')

      vcf_files_annotated = SNPEFF_ONEPERLINE_SNP.out.vcf.join(SNPEFF_ONEPERLINE_INDEL.out.vcf)
      GATK_MERGEVCF_ANNOTATED(vcf_files_annotated, 'SNP_INDEL_filtered_annotated_final')
      
      SNPSIFT_EXTRACTFIELDS(GATK_MERGEVCF_ANNOTATED.out.vcf)

  } else if (params.gen_org=='mouse'){

    // Step 6: Variant Pre-Processing - Part 3
      collecths_metric = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai)
      PICARD_COLLECTHSMETRICS(collecths_metric)
                              
    // Step 7: Variant Calling
      haplotype_caller = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai)
      GATK_HAPLOTYPECALLER(haplotype_caller, 'variant')

    // Step 8: Variant Filtration

      SNPSIFT_ANNOTATE_DBSNP(GATK_HAPLOTYPECALLER.out.vcf, params.dbSNP, params.dbSNP_index, 'intermediate')
      
      GATK_INDEXFEATUREFILE(SNPSIFT_ANNOTATE_DBSNP.out.vcf)

      var_filter = SNPSIFT_ANNOTATE_DBSNP.out.vcf.join(GATK_INDEXFEATUREFILE.out.idx)

      GATK_VARIANTFILTRATION(var_filter, 'BOTH')

    // SNP for final save
      select_var_snp = GATK_VARIANTFILTRATION.out.vcf.join(GATK_VARIANTFILTRATION.out.idx)
      GATK_SELECTVARIANTS_SNP(select_var_snp, 'SNP', 'SNP_filtered_dbsnpID')

    // INDEL for final save
      select_var_indel = GATK_VARIANTFILTRATION.out.vcf.join(GATK_VARIANTFILTRATION.out.idx)
      GATK_SELECTVARIANTS_INDEL(select_var_indel, 'INDEL', 'INDEL_filtered_dbsnpID')

    // Step 9: Post Variant Calling Processing
      SNPEFF(GATK_VARIANTFILTRATION.out.vcf, 'BOTH', 'vcf')

      SNPEFF_ONEPERLINE(SNPEFF.out.vcf, 'BOTH')

      SNPSIFT_EXTRACTFIELDS(SNPEFF_ONEPERLINE.out.vcf)

  }

  agg_stats = JAX_TRIMMER.out.quality_stats.join(PICARD_COLLECTHSMETRICS.out.hsmetrics).join(PICARD_MARKDUPLICATES.out.dedup_metrics)

  // Step 11: Aggregate Stats
  AGGREGATE_STATS(agg_stats)
  
  ch_multiqc_files = Channel.empty()
  ch_multiqc_files = ch_multiqc_files.mix(JAX_TRIMMER.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(ch_GATK_BASERECALIBRATOR_multiqc.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTHSMETRICS.out.hsmetrics.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.dedup_metrics.collect{it[1]}.ifEmpty([]))

  MULTIQC (
      ch_multiqc_files.collect()
  )

}
