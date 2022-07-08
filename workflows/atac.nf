#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/atac.nf'
include {param_log} from '../bin/log/atac.nf'
include {getLibraryId} from '../bin/shared/getLibraryId.nf'
include {CONCATENATE_READS_PE} from '../modules/utility_modules/concatenate_reads_PE'
include {CONCATENATE_READS_SE} from '../modules/utility_modules/concatenate_reads_SE'
include {TRIM_FASTQ} from '../modules/cutadapt/cutadapt_trim_fastq'
include {FASTQC} from '../modules/fastqc/fastqc'
include {ALIGN_TRIMMED_FASTQ} from '../modules/bowtie2/bowtie2_align_trimmed_fastq'
include {SORT as SORT_ALIGN_TRIM;
         SORT as SORT_SHIFTED_BAM;
         SORT as SORT_MARK_DUP_BAM;
         SORT as SORT_LIFTOVER_BAM } from '../modules/samtools/samtools_sort'
include {PICARD_MARKDUPLICATES} from '../modules/picard/picard_markduplicates'
include {REMOVE_DUPLICATE_READS} from '../modules/samtools/samtools_remove_duplicate_reads'
include {CALC_MTDNA_FILTER_CHRM} from '../modules/samtools/samtools_calc_mtdna_filter_chrm'
include {FILTER_REMOVE_MULTI_SHIFT} from '../modules/samtools/samtools_filter_remove_multi_shift'
include {FILTER_REMOVE_MULTI_SIEVE} from '../modules/deeptools/deeptools_filter_remove_multi_sieve'
include {CHAIN_CONVERT_PEAK} from '../modules/g2gtools/g2gtools_chain_convert_peak'
include {CHAIN_EXTRACT_BADREADS} from '../modules/gatk/gatk_chain_extract_badreads'
include {CHAIN_BAD2UNIQ_READS} from '../modules/samtools/samtools_chain_bad2uniq_reads'
include {CHAIN_FILTER_READS} from '../modules/gatk/gatk_chain_filter_reads'
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
include {CALC_PBC_METRICS} from '../modules/bedtools/bedtools_calc_pbc_metrics'
include {LOG_PARSER} from '../modules/python/python_log_parser'

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
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}


// main workflow
workflow ATAC {

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

  // Step 1: Trim_Fastq
  TRIM_FASTQ(read_ch)

  // Step 2: Get fastqc report
  FASTQC(TRIM_FASTQ.out.paired_trimmed_fastq)

  // Step 3: Align trimmed fastq to reference
  ALIGN_TRIMMED_FASTQ(TRIM_FASTQ.out.paired_trimmed_fastq)

  // Step 4: Sort alignment file
  SORT_ALIGN_TRIM(ALIGN_TRIMMED_FASTQ.out.sam, '')

  // Step 5: Flag pcr duplicates
  PICARD_MARKDUPLICATES(SORT_ALIGN_TRIM.out) 

  // Step 6: Remove pcr duplicates
  REMOVE_DUPLICATE_READS(PICARD_MARKDUPLICATES.out.dedup_bam, PICARD_MARKDUPLICATES.out.dedup_bai)

  // Step 7: Calculate %mtDNA and filter mitochondrial reads
  CALC_MTDNA_FILTER_CHRM(REMOVE_DUPLICATE_READS.out.rmDup_bam, REMOVE_DUPLICATE_READS.out.rmDup_bai)

  // Step 8: Filter non-unique and include only 'properly mapped reads' alignments
  FILTER_REMOVE_MULTI_SHIFT(CALC_MTDNA_FILTER_CHRM.out.rmChrM_bam, CALC_MTDNA_FILTER_CHRM.out.rmChrM_bai)

  // Step 9: Run deeptools alignmentSieve 
  FILTER_REMOVE_MULTI_SIEVE(FILTER_REMOVE_MULTI_SHIFT.out[0])

  // Step 10: Re-sort shifted bam
  SORT_SHIFTED_BAM(FILTER_REMOVE_MULTI_SIEVE.out[0], '')


  // If Mouse
  if (params.gen_org=='mouse'){
    // Step 11: Convert peak coordinates
    CHAIN_CONVERT_PEAK(SORT_SHIFTED_BAM.out[0])

    // Step 12: Sort bam by coordinates
    SORT_LIFTOVER_BAM(CHAIN_CONVERT_PEAK.out[0], '')

    // Step 13: Extract a list of 'bad reads'
    CHAIN_EXTRACT_BADREADS(SORT_LIFTOVER_BAM.out[0])

    // Step 14: Remove 'bad reads' from bam file
    CHAIN_BAD2UNIQ_READS(CHAIN_EXTRACT_BADREADS.out.bad_reads)

    // Step 15: Filter list to unique names
    CHAIN_FILTER_READS(SORT_ALIGN_TRIM.out[0], CHAIN_BAD2UNIQ_READS.out.uniq_reads)

    // Step 16: Sort fixmate bam and filter mitochondrial reads
    CHAIN_SORT_FIXMATE_BAM(CHAIN_FILTER_READS.out[0])

    // Step 17: Reference strain samples, filter mitochondrial, unplaced/unlocalized reads and reindex
    NON_CHAIN_REINDEX(SORT_SHIFTED_BAM.out[0])

    // Step 18 : Mix chain and non-chain
    data_ch = CHAIN_SORT_FIXMATE_BAM.out[0].mix(NON_CHAIN_REINDEX.out[0])

  }
  else if (params.gen_org=='human'){
    data_ch = SORT_SHIFTED_BAM.out[0]
  } 

  // Step 19: Peak calling
  PEAK_CALLING(data_ch) 

  if (params.gen_org=='mouse'){
    // Step 20: deeptools bamCoverage bigwig
    BAM_COVERAGE_BIGWIG(data_ch)

  } 

  // Step 21: Fraction of reads in peaks (FRiP)
  FRIP_READS_IN_PEAKS(data_ch, PEAK_CALLING.out.np) 

  // Step 22: Final calculate (FRiP)
  FINAL_CALC_FRIP(data_ch, FRIP_READS_IN_PEAKS.out[0]) 

  // Step 23: Get coverage in each peak
  PEAK_COVERAGE(PEAK_CALLING.out.np) 

  // Step 24: Feature counts
  FEATURE_COUNTS(data_ch, PEAK_COVERAGE.out) 

  if (params.gen_org=='mouse'){
    // Step 25: Feature count to bed
    FEATURE_COUNT2BED(FEATURE_COUNTS.out) 
  }

  // Step 26: Get the fragment length count for quality checks
  QUALITY_CHECKS(FILTER_REMOVE_MULTI_SHIFT.out.srf_bam) 

  // Step 27: Fragment length plot
  FRAG_LEN_PLOT(QUALITY_CHECKS.out) 

  // Step 28: Sort markduplicates bam by read names 
  SORT_MARK_DUP_BAM(PICARD_MARKDUPLICATES.out.dedup_bam, '-n ') 

  // Step 29: Calculating PBC Metrics
  CALC_PBC_METRICS(SORT_MARK_DUP_BAM.out[0]) 

  // Step 30: Log Parser
  LOG_PARSER(TRIM_FASTQ.out.cutadapt_log, ALIGN_TRIMMED_FASTQ.out.bowtie_log, PICARD_MARKDUPLICATES.out.dedup_metrics, CALC_MTDNA_FILTER_CHRM.out.mtdna_log, CALC_PBC_METRICS.out, FINAL_CALC_FRIP.out)

}

