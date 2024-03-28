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
include {CUTADAPT} from "${projectDir}/modules/cutadapt/cutadapt"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {SAMTOOLS_SORT as SAMTOOLS_SORT_PRIMERCLIP;
         SAMTOOLS_SORT as SAMTOOLS_SORT_CALLING} from "${projectDir}/modules/samtools/samtools_sort"
include {PRIMERCLIP} from "${projectDir}/modules/primerclip/primerclip"
include {TARGET_COVERAGE_METRICS} from "${projectDir}/modules/bedtools/bedtools_amplicon_metrics"
include {PICARD_COLLECTTARGETPCRMETRICS} from "${projectDir}/modules/picard/picard_collecttargetpcrmetrics"
include {GATK_BASERECALIBRATOR} from "${projectDir}/modules/gatk/gatk_baserecalibrator"
include {GATK_APPLYBQSR} from "${projectDir}/modules/gatk/gatk_applybqsr"
include {GATK_HAPLOTYPECALLER} from "${projectDir}/modules/gatk/gatk_haplotypecaller"
include {SNPSIFT_ANNOTATE} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {GENERATE_FINGERPRINT_REPORT} from "${projectDir}/modules/python/python_generate_fingerprint_report"
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

  CUTADAPT(read_ch)

  FASTQC(CUTADAPT.out.paired_trimmed_fastq)

  // Step 2: Get Read Group Information
  READ_GROUPS(CUTADAPT.out.paired_trimmed_fastq, "gatk")

  // Step 3: BWA-MEM Alignment
  bwa_mem_mapping = CUTADAPT.out.paired_trimmed_fastq.join(READ_GROUPS.out.read_groups)
                    .map{it -> [it[0], it[1], 'aln', it[2]]}
                      
  BWA_MEM(bwa_mem_mapping)

  SAMTOOLS_SORT_PRIMERCLIP(BWA_MEM.out.sam, '-O sam -n', 'sam')

  PRIMERCLIP(SAMTOOLS_SORT_PRIMERCLIP.out.sorted_file)

  SAMTOOLS_SORT_CALLING(PRIMERCLIP.out.sam, '-O bam', 'bam')

  /*
  Important: While the use of the Picard tool, MarkDuplicates, is a common quality control step to identify
  low-complexity libraries, MarkDuplicates cannot be used on data derived from PCR-based target enrichment
  methods such as the xGen Amplicon Panels. Since these targeted panels contain high numbers of identical
  library fragments (particularly regarding alignment start position), MarkDuplicates cannot appropriately 
  analyze Amplicon libraries.
  https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/application-note/primerclip-a-tool-for-trimming-primer-sequences-application-note.pdf?sfvrsn=cf83e107_14
  */

  GATK_BASERECALIBRATOR(SAMTOOLS_SORT_CALLING.out.sorted_file)
  
  GATK_APPLYBQSR(SAMTOOLS_SORT_CALLING.out.sorted_file.join(GATK_BASERECALIBRATOR.out.table))

  PICARD_COLLECTTARGETPCRMETRICS(GATK_APPLYBQSR.out.bam)

  TARGET_COVERAGE_METRICS(GATK_APPLYBQSR.out.bam)

  GATK_HAPLOTYPECALLER(GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai), '')

  SNPSIFT_ANNOTATE(GATK_HAPLOTYPECALLER.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')

  GENERATE_FINGERPRINT_REPORT(SNPSIFT_ANNOTATE.out.vcf)

  // MultiQC
  ch_multiqc_files = Channel.empty()
  ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.cutadapt_log.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(GATK_BASERECALIBRATOR.out.table.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTTARGETPCRMETRICS.out.txt.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PRIMERCLIP.out.log.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(TARGET_COVERAGE_METRICS.out.qc_metrics.collect{it[1]}.ifEmpty([]))

  MULTIQC (
      ch_multiqc_files.collect()
  )
}