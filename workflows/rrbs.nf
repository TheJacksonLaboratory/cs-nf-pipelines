nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/rrbs'
include {param_log} from '../bin/log/rrbs'
include {getLibraryId} from '../bin/shared/getLibraryId.nf'
include {CONCATENATE_READS_PE} from '../modules/utility_modules/concatenate_reads_PE'
include {CONCATENATE_READS_SE} from '../modules/utility_modules/concatenate_reads_SE'
include {FASTQC} from '../modules/fastqc/fastqc'
include {TRIM_GALORE} from '../modules/trim_galore/trim_galore'
include {BISMARK_ALIGNMENT} from '../modules/bismark/bismark_alignment'
include {BISMARK_DEDUPLICATION} from '../modules/bismark/bismark_deduplication'
include {BISMARK_METHYLATION_EXTRACTION} from '../modules/bismark/bismark_methylation_extraction'
include {MULTIQC} from '../modules/multiqc/multiqc'

// help if needed
if (params.help){
    help()
    exit 0
}

// log paramiter info
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
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

// main workflow
workflow RRBS {

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

  FASTQC(read_ch)

  TRIM_GALORE(read_ch)

  BISMARK_ALIGNMENT(TRIM_GALORE.out.trimmed_fastq)

  BISMARK_DEDUPLICATION(BISMARK_ALIGNMENT.out.bam)

  BISMARK_METHYLATION_EXTRACTION(BISMARK_DEDUPLICATION.out.dedup_bam)


  ch_multiqc_files = Channel.empty()
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(TRIM_GALORE.out.trim_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(TRIM_GALORE.out.trimmed_fastqc.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(BISMARK_ALIGNMENT.out.report.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(BISMARK_DEDUPLICATION.out.dedup_report.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(BISMARK_METHYLATION_EXTRACTION.out.extractor_reports.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(BISMARK_METHYLATION_EXTRACTION.out.extractor_mbias.collect{it[1]}.ifEmpty([]))

  MULTIQC (
      ch_multiqc_files.collect()
  )

}
