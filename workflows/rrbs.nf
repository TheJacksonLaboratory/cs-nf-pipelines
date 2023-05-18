nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/rrbs"
include {param_log} from "${projectDir}/bin/log/rrbs"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv.nf"
include {FILE_DOWNLOAD} from "${projectDir}/subworkflows/aria_download_parse"
include {CONCATENATE_LOCAL_FILES} from "${projectDir}/subworkflows/concatenate_local_files"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {TRIM_GALORE} from "${projectDir}/modules/trim_galore/trim_galore"
include {BISMARK_ALIGNMENT} from "${projectDir}/modules/bismark/bismark_alignment"
include {SAMTOOLS_SORT} from "${projectDir}/modules/samtools/samtools_sort"
include {SAMTOOLS_INDEX} from "${projectDir}/modules/samtools/samtools_index"
include {BISMARK_DEDUPLICATION} from "${projectDir}/modules/bismark/bismark_deduplication"
include {BISMARK_METHYLATION_EXTRACTION} from "${projectDir}/modules/bismark/bismark_methylation_extraction"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

// help if needed
if (params.help){
    help()
    exit 0
}

// log paramiter info
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

// main workflow
workflow RRBS {
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

  FASTQC(read_ch)
  // Note: fastqc is run prior to trimming, as trim galor outputs fastqc level data.

  TRIM_GALORE(read_ch)

  BISMARK_ALIGNMENT(TRIM_GALORE.out.trimmed_fastq)

  SAMTOOLS_SORT(BISMARK_ALIGNMENT.out.bam, '-O bam', 'bam')
  SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sorted_file)

  ch_BISMARK_DEDUPLICATION_multiqc = Channel.empty()

  if (params.skip_deduplication) {
    BISMARK_METHYLATION_EXTRACTION(BISMARK_ALIGNMENT.out.bam)
  } else {
    BISMARK_DEDUPLICATION(BISMARK_ALIGNMENT.out.bam)
    BISMARK_METHYLATION_EXTRACTION(BISMARK_DEDUPLICATION.out.dedup_bam)
    ch_BISMARK_DEDUPLICATION_multiqc = BISMARK_DEDUPLICATION.out.dedup_report
  }

  ch_multiqc_files = Channel.empty()
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(TRIM_GALORE.out.trim_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(TRIM_GALORE.out.trimmed_fastqc.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(BISMARK_ALIGNMENT.out.report.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(ch_BISMARK_DEDUPLICATION_multiqc.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(BISMARK_METHYLATION_EXTRACTION.out.extractor_reports.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(BISMARK_METHYLATION_EXTRACTION.out.extractor_mbias.collect{it[1]}.ifEmpty([]))

  MULTIQC (
      ch_multiqc_files.collect()
  )

}
