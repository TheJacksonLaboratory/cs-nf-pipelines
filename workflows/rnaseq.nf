#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {READ_GROUPS} from '../modules/read_groups'
include {SUMMARY_STATS} from '../modules/summary_stats'
include {BAMTOOLS_RNASEQ_MOUSE} from '../modules/bamtools'
include {RSEM_ALIGNMENT_EXPRESSION} from '../modules/rsem'
include {QUALITY_STATISTICS} from '../modules/quality_stats'
include {PICARD_ALN_METRICS_A;PICARD_ALN_METRICS_B} from '../modules/picard'
include {GATK_STATS_A;GATK_STATS_B} from '../modules/gatk'

// prepare reads channel *
if (params.read_type == 'PE'){
  read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}_*${params.extension}",checkExists:true )
}
else if (params.read_type == 'SE'){
  read_ch = Channel.fromFilePairs("${params.fq_path}/*${params.extension}",checkExists:true, size:1 )
}

// downstream resources (only load once so do it here)
rsem_ref_files = file("${params.rsem_ref_files}/*")

// main workflow
workflow RNASEQ {
  
  // Step 1: Qual_Stat
  QUALITY_STATISTICS(read_ch)

  // Step 2: RSEM
  RSEM_ALIGNMENT_EXPRESSION(QUALITY_STATISTICS.out.trimmed_fastq, rsem_ref_files)

  //Step 3: Get Read Group Information
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq)

  // Step 4a: Picard Alignment Metrics
  PICARD_ALN_METRICS_A(READ_GROUPS.out.read_groups,
                       RSEM_ALIGNMENT_EXPRESSION.out.genome_sorted_bam)
  
  if ("${params.gen_org}" == 'mouse'){
    // Step 4b: Bamtools
    BAMTOOLS_RNASEQ_MOUSE(PICARD_ALN_METRICS_A.out.reordered_sorted_bam)  
  }

  //HUMAN NEEDS TO BE FLUSHED OUT
  if ("${params.gen_org}" == 'human'){

    // Step 4b: Picard Alignment Metrics
    PICARD_ALN_METRICS_B(PICARD_ALN_METRICS_A.out.reordered_sorted_bam)

    // Step 5: Summary Stats
    SUMMARY_STATS(RSEM_ALIGNMENT_EXPRESSION.out.rsem_stats,
                  QUALITY_STATISTICS.out.quality_stats,
                  PICARD_ALN_METRICS_B.out.picard_metrics)

    // Step 6a: GATK Coverage Stats
    GATK_STATS_A(PICARD_ALN_METRICS_A.out.reordered_sorted_bam,
                 PICARD_ALN_METRICS_A.out.reordered_sorted_bai)

    // Step 6b: GATK Coverage Stats
    GATK_STATS_B(GATK_STATS_A.out.gatk_3,
                 GATK_STATS_A.out.gatk_6)
   
  }
}
