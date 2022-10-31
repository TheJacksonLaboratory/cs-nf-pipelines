#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/atac.nf"
include {param_log} from "${projectDir}/bin/log/atac.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {TRIM_FASTQ} from "${projectDir}/modules/cutadapt/cutadapt_trim_fastq"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {ALIGN_TRIMMED_FASTQ} from "${projectDir}/modules/bowtie2/bowtie2_align_trimmed_fastq"
include {SORT as SORT_ALIGN_TRIM;
         SORT as SORT_SHIFTED_BAM;
         SORT as SORT_MARK_DUP_BAM;
         SORT as SORT_LIFTOVER_BAM } from "${projectDir}/modules/samtools/samtools_sort"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {REMOVE_DUPLICATE_READS} from "${projectDir}/modules/samtools/samtools_remove_duplicate_reads"
include {CALC_MTDNA_FILTER_CHRM} from "${projectDir}/modules/samtools/samtools_calc_mtdna_filter_chrm"
include {FILTER_REMOVE_MULTI_SHIFT} from "${projectDir}/modules/samtools/samtools_filter_remove_multi_shift"
include {FILTER_REMOVE_MULTI_SIEVE} from "${projectDir}/modules/deeptools/deeptools_filter_remove_multi_sieve"
include {CHAIN_CONVERT} from "${projectDir}/modules/g2gtools/g2gtools_chain_convert_peak"
include {CHAIN_EXTRACT_BADREADS} from "${projectDir}/modules/gatk/gatk_chain_extract_badreads"
include {CHAIN_BAD2UNIQ_READS} from "${projectDir}/modules/samtools/samtools_chain_bad2uniq_reads"
include {CHAIN_FILTER_READS} from "${projectDir}/modules/gatk/gatk_chain_filter_reads"
include {CHAIN_SORT_FIXMATE_BAM} from "${projectDir}/modules/samtools/samtools_chain_sort_fixmate_bam"
include {NON_CHAIN_REINDEX} from "${projectDir}/modules/samtools/samtools_non_chain_reindex"
include {PEAK_CALLING} from "${projectDir}/modules/macs2/macs2_peak_calling"
include {BAM_COVERAGE_BIGWIG} from "${projectDir}/modules/deeptools/deeptools_bam_coverage_bigwig"
include {FRIP_READS_IN_PEAKS} from "${projectDir}/modules/bedtools/bedtools_frip_reads_in_peaks"
include {FINAL_CALC_FRIP} from "${projectDir}/modules/samtools/samtools_final_calc_frip"
include {PEAK_COVERAGE} from "${projectDir}/modules/macs2/macs2_peak_coverage"
include {FEATURE_COUNTS} from "${projectDir}/modules/subread/subread_feature_counts"
include {FEATURE_COUNT2BED} from "${projectDir}/modules/bedtools/bedtools_feature_count2bed"
include {QUALITY_CHECKS} from "${projectDir}/modules/samtools/samtools_quality_checks"
include {FRAG_LEN_PLOT} from "${projectDir}/modules/rstudio/rstudio_frag_len_plot"
include {CALC_PBC_METRICS} from "${projectDir}/modules/bedtools/bedtools_calc_pbc_metrics"
include {LOG_PARSER} from "${projectDir}/modules/python/python_log_parser"

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
  remove_duplicates = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai)
  REMOVE_DUPLICATE_READS(remove_duplicates)

  // Step 7: Calculate %mtDNA and filter mitochondrial reads
  mt_filter = REMOVE_DUPLICATE_READS.out.rmDup_bam.join(REMOVE_DUPLICATE_READS.out.rmDup_bai)
  CALC_MTDNA_FILTER_CHRM(mt_filter)

  // Step 8: Filter non-unique and include only 'properly mapped reads' alignments
  filter_reads = CALC_MTDNA_FILTER_CHRM.out.rmChrM_bam.join(CALC_MTDNA_FILTER_CHRM.out.rmChrM_bai)
  FILTER_REMOVE_MULTI_SHIFT(filter_reads)

  // Step 9: Run deeptools alignmentSieve 
  FILTER_REMOVE_MULTI_SIEVE(FILTER_REMOVE_MULTI_SHIFT.out[0])

  // Step 10: Re-sort shifted bam
  SORT_SHIFTED_BAM(FILTER_REMOVE_MULTI_SIEVE.out[0], '')


  // If Mouse
  if (params.gen_org=='mouse'){
    // Step 11: Convert peak coordinates
    //          Step occurs when chain != null || chain != false
    CHAIN_CONVERT(SORT_SHIFTED_BAM.out[0])

    // Step 12: Sort bam by coordinates
    SORT_LIFTOVER_BAM(CHAIN_CONVERT.out[0], '')

    // Step 13: Extract a list of 'bad reads'
    CHAIN_EXTRACT_BADREADS(SORT_LIFTOVER_BAM.out[0])

    // Step 14: Remove 'bad reads' from bam file
    CHAIN_BAD2UNIQ_READS(CHAIN_EXTRACT_BADREADS.out.bad_reads)

    // Step 15: Filter list to unique names
    filter_chain_reads = SORT_LIFTOVER_BAM.out[0].join(CHAIN_BAD2UNIQ_READS.out.uniq_reads)
    CHAIN_FILTER_READS(filter_chain_reads)

    // Step 16: Sort fixmate bam and filter mitochondrial reads
    CHAIN_SORT_FIXMATE_BAM(CHAIN_FILTER_READS.out[0])

    // Step 17: Reference strain samples, filter mitochondrial, unplaced/unlocalized reads and reindex
    //          Step occurs when chain == null || chain == false
    NON_CHAIN_REINDEX(SORT_SHIFTED_BAM.out[0])

    // Step 18 : Mix chain and non-chain

    data_ch = CHAIN_SORT_FIXMATE_BAM.out[0].mix(NON_CHAIN_REINDEX.out[0])
    // NOTE: This step will pass only 1 tuple forward [sample_id, [bam, bam.idx]].
    //       The mix statement is required because step 11-16 will only run when `--chain [file]` is called (controlled via modules). 
    //       Step 17 will only run when `--chain` is not used (controlled via modules). 
    //       A bam file is required in the next step. `mix` ensures that one OR the other output is used. 
    //       When '--gen_org == human' data_ch is set to the tuple output in step 10.  

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
  frip_start = data_ch.join(PEAK_CALLING.out.np)
  FRIP_READS_IN_PEAKS(frip_start) 

  // Step 22: Final calculate (FRiP)
  frip_calc = data_ch.join(FRIP_READS_IN_PEAKS.out[0])
  FINAL_CALC_FRIP(frip_calc) 

  // Step 23: Get coverage in each peak
  PEAK_COVERAGE(PEAK_CALLING.out.np) 

  // Step 24: Feature counts
  feature_counts = data_ch.join(PEAK_COVERAGE.out)
  FEATURE_COUNTS(feature_counts) 

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
  log_agg = TRIM_FASTQ.out.cutadapt_log.join(ALIGN_TRIMMED_FASTQ.out.bowtie_log).join(PICARD_MARKDUPLICATES.out.dedup_metrics).join(CALC_MTDNA_FILTER_CHRM.out.mtdna_log).join(CALC_PBC_METRICS.out).join(FINAL_CALC_FRIP.out)
  LOG_PARSER(log_agg)

}

