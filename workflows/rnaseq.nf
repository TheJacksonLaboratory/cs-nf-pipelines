#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/rnaseq"
include {param_log} from "${projectDir}/bin/log/rnaseq"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {XENOME_CLASSIFY} from "${projectDir}/modules/xenome/xenome"
include {FASTQ_SORT as XENOME_SORT} from "${projectDir}/modules/fastq_sort/fastq-tools_sort"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {RNA_SUMMARY_STATS} from "${projectDir}/modules/utility_modules/aggregate_stats_rna"
include {BAMTOOLS_STATS} from "${projectDir}/modules/bamtools/bamtools_stats"
include {RSEM_ALIGNMENT_EXPRESSION} from "${projectDir}/modules/rsem/rsem_alignment_expression"
include {QUALITY_STATISTICS} from "${projectDir}/modules/utility_modules/quality_stats"
include {PICARD_ADDORREPLACEREADGROUPS} from "${projectDir}/modules/picard/picard_addorreplacereadgroups"
include {PICARD_REORDERSAM} from "${projectDir}/modules/picard/picard_reordersam"
include {PICARD_COLLECTRNASEQMETRICS} from "${projectDir}/modules/picard/picard_collectrnaseqmetrics"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"

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

if (params.pdx && params.gen_org == 'mouse') {
    exit 1, "PDX analysis was specified with `--pdx`. `--gen_org` was set to: ${params.gen_org}. This is an invalid parameter combination. `params.gen_org` must == 'human' for PDX analysis."
}

// downstream resources (only load once so do it here)
if (params.rsem_aligner == "bowtie2") {
  rsem_ref_files = file("${params.rsem_ref_files}/bowtie2/*")
}
else if (params.rsem_aligner == "star") {
  rsem_ref_files = file("${params.rsem_ref_files}/STAR/${params.rsem_star_prefix}/*")
}
else error "${params.rsem_aligner} is not valid, use 'bowtie2' or 'star'"

// main workflow
workflow RNASEQ {

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
  QUALITY_STATISTICS(read_ch)

  // Step 1a: Xenome if PDX data used.
  if (params.pdx){
    // Xenome Classification
    XENOME_CLASSIFY(QUALITY_STATISTICS.out.trimmed_fastq)

    // Xenome Read Sort
    XENOME_SORT(XENOME_CLASSIFY.out.xenome_fastq)

    rsem_input = XENOME_SORT.out.sorted_fastq

  } else {
    
    rsem_input = QUALITY_STATISTICS.out.trimmed_fastq

  }


  // Step 2: RSEM
  RSEM_ALIGNMENT_EXPRESSION(rsem_input, rsem_ref_files)

  //Step 3: Get Read Group Information
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq, "picard")

  // Step 4: Picard Alignment Metrics
  add_replace_groups = READ_GROUPS.out.read_groups.join(RSEM_ALIGNMENT_EXPRESSION.out.bam)
  PICARD_ADDORREPLACEREADGROUPS(add_replace_groups)

  PICARD_REORDERSAM(PICARD_ADDORREPLACEREADGROUPS.out.bam)

  // Step 5: Picard Alignment Metrics
  PICARD_SORTSAM(PICARD_REORDERSAM.out.bam)
  // need to sort out ref_flat and ribo_intervals (may break mouse now)
  PICARD_COLLECTRNASEQMETRICS(PICARD_SORTSAM.out.bam)

  // Step 6: Summary Stats

  agg_stats = RSEM_ALIGNMENT_EXPRESSION.out.rsem_stats.join(QUALITY_STATISTICS.out.quality_stats).join(PICARD_COLLECTRNASEQMETRICS.out.picard_metrics)

  RNA_SUMMARY_STATS(agg_stats)

}
