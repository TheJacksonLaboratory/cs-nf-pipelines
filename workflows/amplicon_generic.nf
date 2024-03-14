#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/amplicon.nf"
include {param_log} from "${projectDir}/bin/log/amplicon.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv.nf"
include {FILE_DOWNLOAD} from "${projectDir}/subworkflows/aria_download_parse"
include {CONCATENATE_LOCAL_FILES} from "${projectDir}/subworkflows/concatenate_local_files"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {FASTP} from "${projectDir}/modules/fastp/fastp"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {SAMTOOLS_INDEX} from "${projectDir}/modules/samtools/samtools_index"
include {ABRA2} from "${projectDir}/modules/abra2/abra2"

include {GATK_BASERECALIBRATOR;
         GATK_BASERECALIBRATOR as GATK_BASERECALIBRATOR_ABRA2} from "${projectDir}/modules/gatk/gatk_baserecalibrator"
include {GATK_APPLYBQSR;
         GATK_APPLYBQSR as GATK_APPLYBQSR_ABRA2} from "${projectDir}/modules/gatk/gatk_applybqsr"

include {TARGET_COVERAGE_METRICS} from "${projectDir}/modules/bedtools/bedtools_amplicon_metrics"
include {PICARD_COLLECTTARGETPCRMETRICS} from "${projectDir}/modules/picard/picard_collecttargetpcrmetrics"

include {GATK_HAPLOTYPECALLER} from "${projectDir}/modules/gatk/gatk_haplotypecaller_amplicon"
include {FREEBAYES} from "${projectDir}/modules/freebayes/freebayes"

include {GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from "${projectDir}/modules/gatk/gatk_selectvariants"

include {GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from "${projectDir}/modules/gatk/gatk_variantfiltration"

include {GATK_INDEXFEATUREFILE} from "${projectDir}/modules/gatk/gatk_indexfeaturefile"
include {GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_FREEBAYES} from "${projectDir}/modules/gatk/gatk_variantfiltration_freebayes"

include {BCFTOOLS_NORM as BCFTOOLS_NORM_HAPLOTYPECALLER;
         BCFTOOLS_NORM as BCFTOOLS_NORM_FREEBAYES} from "${projectDir}/modules/bcftools/bcftools_norm"

include {GATK_MERGEVCF} from "${projectDir}/modules/gatk/gatk_mergevcf"

include {ADD_ALT_AF as ADD_ALT_AF_HAPLOTYPECALLER} from "${projectDir}/modules/python/python_add_AF_haplotypecaller"
include {ADD_ALT_AF as ADD_ALT_AF_FREEBAYES} from "${projectDir}/modules/python/python_add_AF_freebayes"

include {BCFTOOLS_MERGECALLERS} from "${projectDir}/modules/bcftools/bcftools_merge_amplicon"

include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_DBSNP;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_COSMIC} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {SNPEFF} from "${projectDir}/modules/snpeff_snpsift/snpeff_snpeff"
include {SNPEFF_ONEPERLINE} from "${projectDir}/modules/snpeff_snpsift/snpeff_oneperline"
include {SNPSIFT_EXTRACTFIELDS} from "${projectDir}/modules/snpeff_snpsift/snpsift_extractfields"
include {SNPSIFT_DBNSFP} from "${projectDir}/modules/snpeff_snpsift/snpsift_dbnsfp"


include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

if (params.download_data && !params.csv_input) {
    exit 1, "Data download was specified with `--download_data`. However, no input CSV file was specified with `--csv_input`. This is an invalid parameter combination. `--download_data` requires a CSV manifest. See `--help` for information."
}

// prepare reads channel
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

workflow AMPLICON {
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
  
  // Step 00: Concat local Fastq files if required.
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

  FASTP(read_ch)

  FASTQC(FASTP.out.trimmed_fastq)

  // Step 2: Get Read Group Information
  READ_GROUPS(FASTP.out.trimmed_fastq, "gatk")

  // Step 3: BWA-MEM Alignment
  bwa_mem_mapping = FASTP.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
  BWA_MEM(bwa_mem_mapping)

  PICARD_SORTSAM(BWA_MEM.out.sam)
  
  if (params.markduplicates) {
    PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
    base_recal_input = PICARD_MARKDUPLICATES.out.dedup_bam
    abra2_input = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai)
  } else {
    base_recal_input = PICARD_SORTSAM.out.bam
    SAMTOOLS_INDEX(PICARD_SORTSAM.out.bam)
    abra2_input = PICARD_SORTSAM.out.bam.join(SAMTOOLS_INDEX.out.bai)
  }

  /*
  Important: While the use of the Picard tool, MarkDuplicates, is a common quality control step to identify
  low-complexity libraries, MarkDuplicates cannot be used on data derived from PCR-based target enrichment
  methods such as the xGen Amplicon Panels. Since these targeted panels contain high numbers of identical
  library fragments (particularly regarding alignment start position), MarkDuplicates cannot appropriately 
  analyze Amplicon libraries.
  https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/application-note/primerclip-a-tool-for-trimming-primer-sequences-application-note.pdf?sfvrsn=cf83e107_14
  */

  // Freebayes Calling: 
  ABRA2(abra2_input)

  GATK_BASERECALIBRATOR_ABRA2(ABRA2.out.bam)

  apply_abra2_bqsr_input = ABRA2.out.bam.join(GATK_BASERECALIBRATOR_ABRA2.out.table)

  GATK_APPLYBQSR_ABRA2(apply_abra2_bqsr_input)
  
  FREEBAYES(GATK_APPLYBQSR_ABRA2.out.bam.join(GATK_APPLYBQSR_ABRA2.out.bai), '')
  
  GATK_INDEXFEATUREFILE(FREEBAYES.out.vcf)

  // HaplotypeCaller Calling: 
  GATK_BASERECALIBRATOR(base_recal_input)

  if (params.markduplicates) {
    apply_bqsr_input = PICARD_MARKDUPLICATES.out.dedup_bam.join(GATK_BASERECALIBRATOR.out.table)
  } else {
    apply_bqsr_input = PICARD_SORTSAM.out.bam.join(GATK_BASERECALIBRATOR.out.table)
  }

  GATK_APPLYBQSR(apply_bqsr_input)

  PICARD_COLLECTTARGETPCRMETRICS(GATK_APPLYBQSR.out.bam)

  TARGET_COVERAGE_METRICS(GATK_APPLYBQSR.out.bam)

  GATK_HAPLOTYPECALLER(GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai), '')
  

  // Step 8: Variant Filtration

  // Freebayes
  var_filter = FREEBAYES.out.vcf.join(GATK_INDEXFEATUREFILE.out.idx)
  GATK_VARIANTFILTRATION_FREEBAYES(var_filter, 'freebayes_filtered.vcf')

  BCFTOOLS_NORM_FREEBAYES(GATK_VARIANTFILTRATION_FREEBAYES.out.vcf)

  ADD_ALT_AF_FREEBAYES(BCFTOOLS_NORM_FREEBAYES.out.vcf, 'freebayes_altAF_filtered')


  // HaplotypeCaller
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

  vcf_files = GATK_VARIANTFILTRATION_SNP.out.vcf.join(GATK_VARIANTFILTRATION_INDEL.out.vcf)
  GATK_MERGEVCF (vcf_files, 'haplotypecaller_filtered')

  BCFTOOLS_NORM_HAPLOTYPECALLER(GATK_MERGEVCF.out.vcf)

  ADD_ALT_AF_HAPLOTYPECALLER(BCFTOOLS_NORM_HAPLOTYPECALLER.out.vcf, 'haplotypecaller_altAF_filtered')

  vcf_merge_input = ADD_ALT_AF_HAPLOTYPECALLER.out.vcf.join(ADD_ALT_AF_FREEBAYES.out.vcf)

  BCFTOOLS_MERGECALLERS(vcf_merge_input, 'mergedCallers_filtered_unannotated')

  // Step 9: Variant Annotation

  SNPSIFT_ANNOTATE_DBSNP(BCFTOOLS_MERGECALLERS.out.vcf, params.dbSNP, params.dbSNP_index, 'dbSNP')
  SNPSIFT_ANNOTATE_COSMIC(SNPSIFT_ANNOTATE_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmic')
  SNPEFF(SNPSIFT_ANNOTATE_COSMIC.out.vcf, 'BOTH', 'vcf')
  SNPSIFT_DBNSFP(SNPEFF.out.vcf, 'BOTH')
  SNPEFF_ONEPERLINE(SNPSIFT_DBNSFP.out.vcf, 'BOTH')

  SNPSIFT_EXTRACTFIELDS(SNPEFF_ONEPERLINE.out.vcf)

  // MultiQC
  ch_multiqc_files = Channel.empty()
  ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.quality_json.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(GATK_BASERECALIBRATOR.out.table.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTTARGETPCRMETRICS.out.txt.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(TARGET_COVERAGE_METRICS.out.qc_metrics.collect{it[1]}.ifEmpty([]))

  MULTIQC (
      ch_multiqc_files.collect()
  )
}
