#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/rnaseq'
include {param_log} from '../bin/log/rnaseq'
include {READ_GROUPS} from '../modules/read_groups'
include {SUMMARY_STATS} from '../modules/summary_stats'
include {BAMTOOLS_STATS} from '../modules/bamtools'
include {RSEM_ALIGNMENT_EXPRESSION} from '../modules/rsem'
include {QUALITY_STATISTICS} from '../modules/quality_stats'
include {PICARD_ADDORREPLACEREADGROUPS;
         PICARD_REORDERSAM;
         PICARD_COLLECTRNASEQMETRICS;
         PICARD_SORTSAM} from '../modules/picard'
include {GATK_DEPTHOFCOVERAGE as GATK_DEPTHOFCOVERAGE_CTP;
         GATK_DEPTHOFCOVERAGE as GATK_DEPTHOFCOVERAGE_PROBES} from '../modules/gatk'
include {FORMAT_GATK as FORMAT_GATK_CTP;
         FORMAT_GATK as FORMAT_GATK_PROBES;
         COVCALC_GATK as COVCALC_GATK_CTP;
         COVCALC_GATK as COVCALC_GATK_PROBES} from '../bin/rnaseq/gatk_formatter'

// help if needed
if (params.help){
    help()
    exit 0
}

// log paramiter info
param_log()

// prepare reads channel *
if (params.read_type == 'PE'){
  read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
}
else if (params.read_type == 'SE'){
  read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
}

// if channel is empty give error message and exit
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

// downstream resources (only load once so do it here)
rsem_ref_files = file("${params.rsem_ref_files}/*")

// main workflow
workflow RNASEQ {

  // Step 1: Qual_Stat
  QUALITY_STATISTICS(read_ch)

  // Step 2: RSEM
  RSEM_ALIGNMENT_EXPRESSION(QUALITY_STATISTICS.out.trimmed_fastq, rsem_ref_files)

  //Step 3: Get Read Group Information
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq, "picard")

  // Step 4: Picard Alignment Metrics
  PICARD_ADDORREPLACEREADGROUPS(READ_GROUPS.out.read_groups,
                                RSEM_ALIGNMENT_EXPRESSION.out.bam)
  PICARD_REORDERSAM(PICARD_ADDORREPLACEREADGROUPS.out.bam)

  // Step 5: Picard Alignment Metrics
  PICARD_SORTSAM(PICARD_REORDERSAM.out.bam)
  // need to sort out ref_flat and ribo_intervals (may break mouse now)
  PICARD_COLLECTRNASEQMETRICS(PICARD_SORTSAM.out.bam)

  // Step 6: Summary Stats
  SUMMARY_STATS(RSEM_ALIGNMENT_EXPRESSION.out.rsem_stats,
                QUALITY_STATISTICS.out.quality_stats,
                PICARD_COLLECTRNASEQMETRICS.out.picard_metrics)

  /* do we need this?
  if ("${params.gen_org}" == 'mouse'){
     BAMTOOLS_STATS(PICARD_REORDERSAM.out.bam)
  }
  */

  // If gen_org human
  if ("${params.gen_org}" == 'human'){
    // Step 7: GATK Coverage Stats
      // CTP
        GATK_DEPTHOFCOVERAGE_CTP(PICARD_SORTSAM.out.bam, PICARD_SORTSAM.out.bai, params.ctp_genes)
        FORMAT_GATK_CTP(GATK_DEPTHOFCOVERAGE_CTP.out.txt, params.ctp_genes)
        COVCALC_GATK_CTP(FORMAT_GATK_CTP.out.txt, "CTP")

      // PROBES
        GATK_DEPTHOFCOVERAGE_PROBES(PICARD_SORTSAM.out.bam, PICARD_SORTSAM.out.bai, params.probes)
        FORMAT_GATK_PROBES(GATK_DEPTHOFCOVERAGE_PROBES.out.txt, params.probes)
        COVCALC_GATK_PROBES(FORMAT_GATK_CTP.out.txt, "PROBES")

  }
}
