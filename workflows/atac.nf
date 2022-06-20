#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {TRIM_FASTQ} from '../modules/cutadapt/cutadapt_trim_fastq'
include {FASTQC} from '../modules/fastqc/fastqc'
include {ALIGN_TRIMMED_FASTQ} from '../modules/bowtie2/bowtie2_align_trimmed_fastq'
include {SORT_ALIGNMENT} from '../modules/samtools/samtools_sort_alignment'
include {FLAG_PCR_DUPES} from '../modules/gatk/gatk_flag_pcr_dupes'
include {RM_DUPES_READS} from '../modules/samtools/samtools_rm_dupe_reads'
include {CALC_MTDNA_FILTER_CHRM} from '../modules/samtools/samtools_calc_mtdna_filter_chrm'
include {FILTER_RMMULTI_SHIFT} from '../modules/samtools/samtools_filter_rmmulti_shift'
include {FILTER_RMMULTI_SIEVE} from '../modules/deeptools/deeptools_filter_rmmulti_sieve'
include {FILTER_RMMULTI_SORT} from '../modules/samtools/samtools_filter_rmmulti_sort'
include {CHAIN_CONVERT_PEAK_B6} from '../modules/g2gtools/g2gtools_chain_convert_peak_b6'
include {CHAIN_SORT_COORDS} from '../modules/samtools/samtools_chain_sort_coords'
include {CHAIN_EXTRACT_BADREADS_B6} from '../modules/gatk/gatk_chain_extract_badreads_b6'
include {CHAIN_BAD2UNIQ_READS} from '../modules/samtools/samtools_chain_bad2uniq_reads'
include {CHAIN_FILTER_READS_B6} from '../modules/gatk/gatk_chain_filter_reads_b6'
include {CHAIN_SORT_FIXMATE_BAM} from '../modules/samtools/samtools_chain_sort_fixmate_bam'
include {NON_CHAIN_REINDEX} from '../modules/samtools/samtools_non_chain_reindex'
include {PEAK_CALLING} from '../modules/macs2/macs2_peak_calling'
include {BAM_COVERAGE_BIGWIG} from '../modules/deeptools/deeptools_bam_coverage_bigwig'
include {FRIP_READS_IN_PEAKS} from '../modules/bedtools/bedtools_frip_reads_in_peaks'
include {FINAL_CALC_FRIP} from '../modules/samtools/samtools_final_calc_frip'
include {PEAK_COVERAGE} from '../modules/macs2/macs2_peak_coverage'
include {FEATURE_COUNTS} from '../modules/subread/subread_feature_counts'
include {FEATURE_COUNT2BED} from '../modules/bedtools/bedtools_feature_count2bed'
include {QUALITY_CHECKS} from '../modules/samtools/samtools_quality_checks'
include {FRAG_LEN_PLOT} from '../modules/rstudio/rstudio_frag_len_plot'
include {LIBRARY_COMPLEXITY} from '../modules/samtools/samtools_library_complexity'
include {CALC_PBC_METRICS} from '../modules/bedtools/bedtools_calc_pbc_metrics'
include {LOG_PARSER} from '../modules/python/python_log_parser'


// prepare reads channel
read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )

// if channel is empty give error message and exit
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

// main workflow
workflow ATAC {

  // Step 1: Trim_Fastq
  TRIM_FASTQ(read_ch)

  // Step 2: Get fastqc report
  FASTQC(TRIM_FASTQ.out.paired_trimmed_fastq)

  // Step 3: Align trimmed fastq to reference
  ALIGN_TRIMMED_FASTQ(TRIM_FASTQ.out.paired_trimmed_fastq)

  // Step 4: Sort Alignment file
  SORT_ALIGNMENT(ALIGN_TRIMMED_FASTQ.out.sam)

  // Step 5: Flag PCR Duplicates
  FLAG_PCR_DUPES(SORT_ALIGNMENT.out.bam, SORT_ALIGNMENT.out.bai)

  // Step 6: Remove PCR Duplicates
  RM_DUPES_READS(FLAG_PCR_DUPES.out.marked_bam, FLAG_PCR_DUPES.out.marked_bai)

  // Step 7: Calculate and filter
  CALC_MTDNA_FILTER_CHRM(RM_DUPES_READS.out.rmDup_bam, RM_DUPES_READS.out.rmDup_bai)

  // Step 8: Filter Non-Unique and Include Only 'properly mapped reads' Alignments
  FILTER_RMMULTI_SHIFT(CALC_MTDNA_FILTER_CHRM.out.rmChrM_bam, CALC_MTDNA_FILTER_CHRM.out.rmChrM_bai)

  // Step 9: Run deeptools alignmentSieve 
  FILTER_RMMULTI_SIEVE(FILTER_RMMULTI_SHIFT.out[0])

  // Step 10: Re-sorting shifted bam
  FILTER_RMMULTI_SORT(FILTER_RMMULTI_SIEVE.out[0])


  // If Mouse
  if (params.gen_org=='mouse'){
    // Step 11: Converting Peak Coordinates to B6
    CHAIN_CONVERT_PEAK_B6(FILTER_RMMULTI_SORT.out[0])

    // Step 12: Sorting bam by Coordinates
    CHAIN_SORT_COORDS(CHAIN_CONVERT_PEAK_B6.out[0])

    // Step 13: Extracting a list of 'bad reads'
    CHAIN_EXTRACT_BADREADS_B6(CHAIN_SORT_COORDS.out[0])

    // Step 14: Removing 'bad reads' from bam file
    CHAIN_BAD2UNIQ_READS(CHAIN_EXTRACT_BADREADS_B6.out.bad_reads)

    // Step 15: Filtering list to unique names
    CHAIN_FILTER_READS_B6(SORT_ALIGNMENT.out[0], CHAIN_BAD2UNIQ_READS.out.uniq_reads)

    // Step 16: Performing name sort the bam
    CHAIN_SORT_FIXMATE_BAM(CHAIN_FILTER_READS_B6.out[0])

    // Step 17: Non chain reindex
    NON_CHAIN_REINDEX(CHAIN_SORT_FIXMATE_BAM.out[0])

    // Step 18 : Mixing chain and non-chain
    //data_ch = CHAIN_SORT_FIXMATE_BAM.out[0].mix(NON_CHAIN_REINDEX.out[0])
    data_ch = CHAIN_SORT_FIXMATE_BAM.out.mix(NON_CHAIN_REINDEX.out)

  }
  else if (params.gen_org=='human'){
    data_ch = FILTER_RMMULTI_SORT.out[0]
  } 

  // Step 19: Performing Peak Calling
  PEAK_CALLING(data_ch) 

  if (params.gen_org=='mouse'){
    // Step 20: Running deeptools bamCoverage bigwig
    BAM_COVERAGE_BIGWIG(data_ch)

  } 

  // Step 21: Running Fraction of reads in peaks (FRiP)
  FRIP_READS_IN_PEAKS(data_ch, PEAK_CALLING.out.np) 

  // Step 22: Final Calculate (FRiP)
  FINAL_CALC_FRIP(data_ch, FRIP_READS_IN_PEAKS.out[0]) 

  // Step 23: Get coverage in each peak
  PEAK_COVERAGE(PEAK_CALLING.out.np) 

  // Step 24: Performing Feature Counts
  FEATURE_COUNTS(data_ch, PEAK_COVERAGE.out) 

  if (params.gen_org=='mouse'){
    // Step 25: Performing Feature Count to Bed
    FEATURE_COUNT2BED(FEATURE_COUNTS.out) 
  }

  // Step 26: Performing Feature Count to Bed
  QUALITY_CHECKS(FILTER_RMMULTI_SHIFT.out.srf_bam) 

  // Step 27: Performing Fragment Length Plot
  FRAG_LEN_PLOT(QUALITY_CHECKS.out) 

  // Step 28: Performing Library Complexity 
  LIBRARY_COMPLEXITY(FLAG_PCR_DUPES.out.marked_bam, FLAG_PCR_DUPES.out.marked_bai) 

  // Step 29: Calculating PBC Metrics
  CALC_PBC_METRICS(LIBRARY_COMPLEXITY.out) 

  // Step 30: Log Parser
  LOG_PARSER(TRIM_FASTQ.out.cutadapt_log, ALIGN_TRIMMED_FASTQ.out.bowtie_log, FLAG_PCR_DUPES.out.srt_metr_log, CALC_MTDNA_FILTER_CHRM.out.mtdna_log, CALC_PBC_METRICS.out, FINAL_CALC_FRIP.out)

}

