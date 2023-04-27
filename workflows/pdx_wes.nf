#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/pdx_wes.nf"
include {param_log} from "${projectDir}/bin/log/pdx_wes.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv.nf"
include {FILE_DOWNLOAD} from "${projectDir}/subworkflows/aria_download_parse"
include {CONCATENATE_LOCAL_FILES} from "${projectDir}/subworkflows/concatenate_local_files"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {JAX_TRIMMER} from "${projectDir}/modules/utility_modules/jax_trimmer"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {XENOME_CLASSIFY} from "${projectDir}/modules/xenome/xenome"
include {FASTQ_PAIR} from "${projectDir}/modules/fastq-tools/fastq-pair"
include {FASTQ_SORT as FASTQ_SORT_HUMAN;
         FASTQ_SORT as FASTQ_SORT_MOUSE} from "${projectDir}/modules/fastq-tools/fastq-sort"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {GATK_BASERECALIBRATOR} from "${projectDir}/modules/gatk/gatk_baserecalibrator"
include {GATK_APPLYBQSR} from "${projectDir}/modules/gatk/gatk_applybqsr"
include {GATK_GETSAMPLENAME} from "${projectDir}/modules/gatk/gatk_getsamplename_noMeta"
include {GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from "${projectDir}/modules/gatk/gatk_variantfiltration_mutect2"
include {GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from "${projectDir}/modules/gatk/gatk_selectvariants"
include {GATK_MUTECT2} from "${projectDir}/modules/gatk/gatk_mutect2_tumorOnly"
include {GATK_FILTERMUECTCALLS} from "${projectDir}/modules/gatk/gatk_filtermutectcalls_tumorOnly"
include {MSISENSOR2_MSI} from "${projectDir}/modules/msisensor2/msisensor2_tumorOnly"
include {GATK_MERGEVCF as GATK_MERGEVCF_UNANNOTATED;
         GATK_MERGEVCF as GATK_MERGEVCF_ANNOTATED} from "${projectDir}/modules/gatk/gatk_mergevcf"
include {COSMIC_ANNOTATION as COSMIC_ANNOTATION_SNP;
         COSMIC_ANNOTATION as COSMIC_ANNOTATION_INDEL} from "${projectDir}/modules/cosmic/cosmic_annotation"
include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_DBSNP;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_DBSNP} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {SNPEFF as SNPEFF_SNP;
         SNPEFF as SNPEFF_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_snpeff"
include {SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_SNP;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_oneperline"
include {SNPSIFT_EXTRACTFIELDS} from "${projectDir}/modules/snpeff_snpsift/snpsift_extractfields"
include {SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpsift_dbnsfp"
include {PICARD_COLLECTHSMETRICS} from "${projectDir}/modules/picard/picard_collecthsmetrics"
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

if (params.download_data && !params.csv_input) {
    exit 1, "Data download was specified with `--download_data`. However, no input CSV file was specified with `--csv_input`. This is an invalid parameter combination. `--download_data` requires a CSV manifest. See `--help` for information."
}

if (params.gen_org == 'mouse') {
    exit 1, "PDX workflow was called; however, `--gen_org` was set to: ${params.gen_org}. This is an invalid parameter combination. `--gen_org` must == 'human' for PDX analysis."
}

if (params.csv_input) {

    ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))
    
    if (params.read_type == 'PE'){
        ch_input_sample.map{it -> [it[0], [it[2], it[3]]]}.set{read_ch}
        ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
    } else if (params.read_type == 'SE') {
        ch_input_sample.map{it -> [it[0], it[2]]}.set{read_ch}
        ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

} else if (params.concat_lanes){
  
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
    // if channel is empty give error message and exit
    read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern} and file extension: ${params.extension}"}

} else {
  
  if (params.read_type == 'PE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
  }
    // if channel is empty give error message and exit
    read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern} and file extension: ${params.extension}"}

}

workflow PDX_WES {

    // Step 0: Download data and concat Fastq files if needed. 
    if (params.download_data){
        FILE_DOWNLOAD(ch_input_sample)

        FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[2]]}.set{read_ch}
        FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

    // Step 00: Concat local Fastq files from CSV input if required.
    if (!params.download_data && params.csv_input){
        CONCATENATE_LOCAL_FILES(ch_input_sample)
        CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[2]]}.set{read_ch}
        CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

    // Step 00: Concat local Fastq files from directory if required.
    if (params.concat_lanes && !params.csv_input){
        if (params.read_type == 'PE'){
            CONCATENATE_READS_PE(read_ch)
            read_ch = CONCATENATE_READS_PE.out.concat_fastq
        } else if (params.read_type == 'SE'){
            CONCATENATE_READS_SE(read_ch)
            read_ch = CONCATENATE_READS_SE.out.concat_fastq
        }
    }

    // ** MAIN workflow starts: 

    // Step 1: Qual_Stat
    JAX_TRIMMER(read_ch)

    if (params.read_type == 'PE') {
      FASTQ_PAIR(JAX_TRIMMER.out.trimmed_fastq)
      xenome_input = FASTQ_PAIR.out.paired_fastq
    } else {
      xenome_input = JAX_TRIMMER.out.trimmed_fastq
    }

    FASTQC(JAX_TRIMMER.out.trimmed_fastq)

    // Step 2: Xenome classify and sort. 
    XENOME_CLASSIFY(xenome_input)

    // Xenome Read Sort
    FASTQ_SORT_HUMAN(XENOME_CLASSIFY.out.xenome_fastq, 'human')
    FASTQ_SORT_MOUSE(XENOME_CLASSIFY.out.xenome_mouse_fastq, 'mouse')

    // Step 3: Get Read Group Information
    READ_GROUPS(JAX_TRIMMER.out.trimmed_fastq, "gatk")

    // Step 4: BWA-MEM Alignment
    bwa_mem_mapping = FASTQ_SORT_HUMAN.out.sorted_fastq.join(READ_GROUPS.out.read_groups)

    BWA_MEM(bwa_mem_mapping)

    // Step 5: Variant Preprocessing - Part 1
    PICARD_SORTSAM(BWA_MEM.out.sam)
    PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)

    // Step 6: Variant Pre-Processing - Part 2
    GATK_BASERECALIBRATOR(PICARD_MARKDUPLICATES.out.dedup_bam)

    apply_bqsr = PICARD_MARKDUPLICATES.out.dedup_bam.join(GATK_BASERECALIBRATOR.out.table)
    GATK_APPLYBQSR(apply_bqsr)

    // Step 7: Variant Pre-Processing - Part 3
    collect_metrics = GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai)
    PICARD_COLLECTHSMETRICS(collect_metrics)

    // Step 8: MSI
    MSISENSOR2_MSI(GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai))

    // Step 9: Get sample names
    GATK_GETSAMPLENAME(collect_metrics)

    // ** Variant Calling
    mutect2_caller_input = GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai).join(GATK_GETSAMPLENAME.out.sample_name)
    
    // Step 10: Mutect2
    GATK_MUTECT2(mutect2_caller_input)
    GATK_FILTERMUECTCALLS(GATK_MUTECT2.out.vcf_tbi_stats)

    // Step 8: Variant Filtration
    // SNP
    GATK_SELECTVARIANTS_SNP(GATK_FILTERMUECTCALLS.out.mutect2_vcf_tbi, 'SNP', 'selected_SNP')

    var_filter_snp = GATK_SELECTVARIANTS_SNP.out.vcf.join(GATK_SELECTVARIANTS_SNP.out.idx)
    GATK_VARIANTFILTRATION_SNP(var_filter_snp, 'SNP')

    // INDEL
    GATK_SELECTVARIANTS_INDEL(GATK_FILTERMUECTCALLS.out.mutect2_vcf_tbi, 'INDEL', 'selected_INDEL')

    var_filter_indel = GATK_SELECTVARIANTS_INDEL.out.vcf.join(GATK_SELECTVARIANTS_INDEL.out.idx)
    GATK_VARIANTFILTRATION_INDEL(var_filter_indel, 'INDEL')

    // Step 9: Post Variant Calling Processing - Part 1
    // 
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

    agg_stats = JAX_TRIMMER.out.quality_stats.join(PICARD_COLLECTHSMETRICS.out.hsmetrics).join(PICARD_MARKDUPLICATES.out.dedup_metrics)

    // Step 11: Aggregate Stats
    AGGREGATE_STATS(agg_stats)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(JAX_TRIMMER.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GATK_BASERECALIBRATOR.out.table.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTHSMETRICS.out.hsmetrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.dedup_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(XENOME_CLASSIFY.out.xenome_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GATK_FILTERMUECTCALLS.out.stats.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )

}
