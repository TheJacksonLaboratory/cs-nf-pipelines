#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/chipseq.nf"
include {param_log} from "${projectDir}/bin/log/chipseq.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CHECK_DESIGN} from "${projectDir}/modules/utility_modules/chipseq_check_design"
include {SAMTOOLS_FAIDX} from "${projectDir}/modules/samtools/samtools_faidx"
include {MAKE_GENOME_FILTER} from "${projectDir}/modules/utility_modules/chipseq_make_genome_filter"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {TRIM_GALORE} from "${projectDir}/modules/trim_galore/trim_galore"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {SAMTOOLS_FILTER} from "${projectDir}/modules/samtools/samtools_filter"
include {SAMTOOLS_SORT;
         SAMTOOLS_SORT as PAIR_SORT;
         SAMTOOLS_SORT as NAME_SORT} from "${projectDir}/modules/samtools/samtools_sort"
include {SAMTOOLS_INDEX} from "${projectDir}/modules/samtools/samtools_index"
include {SAMTOOLS_STATS;
         SAMTOOLS_STATS as SAMTOOLS_STATS_MD;
         SAMTOOLS_STATS as SAMTOOLS_STATS_FILTERED;
         SAMTOOLS_STATS as SAMTOOLS_STATS_BF} from "${projectDir}/modules/samtools/samtools_stats"
include {PICARD_MERGESAMFILES} from "${projectDir}/modules/picard/picard_mergesamfiles"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {SAMTOOLS_MERGEBAM_FILTER} from "${projectDir}/modules/samtools/samtools_mergebam_filter"
include {BAMTOOLS_FILTER} from "${projectDir}/modules/bamtools/bamtools_filter"
include {BAMPE_RM_ORPHAN} from "${projectDir}/modules/utility_modules/chipseq_bampe_rm_orphan"
include {PRESEQ} from "${projectDir}/modules/preseq/preseq"
include {PICARD_COLLECTMULTIPLEMETRICS} from "${projectDir}/modules/picard/picard_collectmultiplemetrics"
include {BEDTOOLS_GENOMECOV} from "${projectDir}/modules/bedtools/bedtools_genomecov"
include {UCSC_BEDGRAPHTOBIGWIG} from "${projectDir}/modules/ucsc/ucsc_bedgraphtobigwig"
include {DEEPTOOLS_COMPUTEMATRIX} from "${projectDir}/modules/deeptools/deeptools_computematrix"
include {DEEPTOOLS_PLOTPROFILE} from "${projectDir}/modules/deeptools/deeptools_plotprofile"
include {DEEPTOOLS_PLOTHEATMAP} from "${projectDir}/modules/deeptools/deeptools_plotheatmap"
include {PHANTOMPEAKQUALTOOLS} from "${projectDir}/modules/r/phantompeakqualtools"
include {MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS} from "${projectDir}/modules/r/multiqc_custom_phantompeakqualtools"
include {DEEPTOOLS_PLOTFINGERPRINT} from "${projectDir}/modules/deeptools/deeptools_plotfingerprint"
include {PEAK_CALLING_CHIPSEQ} from "${projectDir}/modules/macs2/macs2_peak_calling_chipseq"
include {FRIP_SCORE} from "${projectDir}/modules/utility_modules/frip_score"
include {HOMER_ANNOTATEPEAKS;
         HOMER_ANNOTATEPEAKS as CONSENSUS_PEAKS_ANNOTATE} from "${projectDir}/modules/homer/homer_annotatepeaks"
include {PLOT_MACS2_QC} from "${projectDir}/modules/macs2/plot_macs2_qc"
include {PLOT_HOMER_ANNOTATEPEAKS} from "${projectDir}/modules/homer/plot_homer_annotatepeaks"
include {MACS2_CONSENSUS} from "${projectDir}/modules/macs2/macs2_consensus"
include {ANNOTATE_BOOLEAN_PEAKS} from "${projectDir}/modules/homer/annotate_boolean_peaks"
include {SUBREAD_FEATURECOUNTS} from "${projectDir}/modules/subread/subread_feature_counts_chipseq"
include {DESEQ2_QC} from "${projectDir}/modules/utility_modules/deseq2_qc"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

// help if needed
if (params.help){
    help()
    exit 0
}

ANSI_RED = "\u001B[31m";
ANSI_RESET = "\u001B[0m";

// log params
param_log()

// main workflow
workflow CHIPSEQ {

  if (params.input)     { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Samples design file not specified!' }

  // Step 1: CHECK_DESIGN
  CHECK_DESIGN(ch_input)

  /*
  * Create channels for input fastq files
  */

  if (params.read_type == 'SE'){
    read_ch = CHECK_DESIGN.out.sample_reads
                  .splitCsv(header:true, sep:',')
                  .map { row -> if (row.fastq_2){
                    System.err.println(ANSI_RED + "---------------------------------------------------------------------------------" + ANSI_RESET)
                    System.err.println(ANSI_RED + "ERROR: Param: `--read_type = SE` was specified, but the csv input contains `fastq_2`. Adjust input csv or parameter and restart." + ANSI_RESET)
                    System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                    System.err.println(ANSI_RED + "---------------------------------------------------------------------------------" + ANSI_RESET)
                    System.exit(1)
                  } else if (!(row.fastq_1)) {
                    System.err.println(ANSI_RED + "---------------------------------------------------------------------------------" + ANSI_RESET)
                    System.err.println(ANSI_RED + "ERROR: Missing field in csv file header. Param: `--read_type = SE` was specified. The csv input file must have field: 'fastq_1'." + ANSI_RESET)
                    System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                    System.err.println(ANSI_RED + "---------------------------------------------------------------------------------" + ANSI_RESET)
                    System.exit(1)
                  }
                  row}
                  .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true) ] ] }
  } else {
      read_ch = CHECK_DESIGN.out.sample_reads
              .splitCsv(header:true, sep:',')
              .map { row -> if (!(row.fastq_1) | !(row.fastq_2)){
                System.err.println(ANSI_RED + "---------------------------------------------------------------------------------" + ANSI_RESET)
                System.err.println(ANSI_RED + "ERROR: Missing field in csv file header. Param: `--read_type = PE` was specified. The csv input file must have fields: 'fastq_1', fastq_2'." + ANSI_RESET)
                System.err.println(ANSI_RED + "Exiting now." + ANSI_RESET)
                System.err.println(ANSI_RED + "---------------------------------------------------------------------------------" + ANSI_RESET)
                System.exit(1)
              }
              row}
            .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
  }

  /*
   * Create a channel with [sample_id, control id, antibody, replicatesExist, multipleGroups]
   */
  control_ch = CHECK_DESIGN.out.study_design
      .splitCsv(header:true, sep:',')
      .map { row -> [ row.sample_id, row.control_id, row.antibody, row.replicatesExist.toBoolean(), row.multipleGroups.toBoolean() ] }


  // Header files for MultiQC
  ch_spp_nsc_header           = file("${projectDir}/bin/shared/multiqc/chipseq/spp_nsc_header.txt", checkIfExists: true)
  ch_spp_rsc_header           = file("${projectDir}/bin/shared/multiqc/chipseq/spp_rsc_header.txt", checkIfExists: true)
  ch_spp_correlation_header   = file("${projectDir}/bin/shared/multiqc/chipseq/spp_correlation_header.txt", checkIfExists: true)
  ch_peak_count_header        = file("${projectDir}/bin/shared/multiqc/chipseq/peak_count_header.txt", checkIfExists: true)
  ch_frip_score_header        = file("${projectDir}/bin/shared/multiqc/chipseq/frip_score_header.txt", checkIfExists: true)
  ch_peak_annotation_header   = file("${projectDir}/bin/shared/multiqc/chipseq/peak_annotation_header.txt", checkIfExists: true)
  ch_deseq2_pca_header        = file("${projectDir}/bin/shared/multiqc/chipseq/deseq2_pca_header.txt", checkIfExists: true)
  ch_deseq2_clustering_header = file("${projectDir}/bin/shared/multiqc/chipseq/deseq2_clustering_header.txt", checkIfExists: true)

  // Reference genome
  ch_fasta = file(params.ref_fa, checkIfExists: true)
  ch_gtf   = file(params.gtf, checkIfExists: true)

  // genes.bed
  if (params.gene_bed)  { ch_gene_bed = file(params.gene_bed, checkIfExists: true) }

  // Step 2: Make genome filter
  faidx_input = ['primary_ref_fasta', ch_fasta]
  SAMTOOLS_FAIDX(faidx_input)
  MAKE_GENOME_FILTER(SAMTOOLS_FAIDX.out.fai, params.blacklist)

  // Step 3: Fastqc
  FASTQC(read_ch)
  
  // Step 4: Trim Galore
  TRIM_GALORE(read_ch)

  // Step 5: Get Read Group Information
  READ_GROUPS(TRIM_GALORE.out.trimmed_fastq, "gatk")

  // Step 6: BWA-MEM
  bwa_mem_mapping = TRIM_GALORE.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
                    .map{it -> [it[0], it[1], 'aln', it[2]]}

  BWA_MEM(bwa_mem_mapping)

  // Step 7: Samtools Removing Unmapped
  SAMTOOLS_FILTER(BWA_MEM.out, '-F 0x0100')

  // Step 8: Samtools Sort
  SAMTOOLS_SORT(SAMTOOLS_FILTER.out.bam, '-O bam', 'bam')

  // Step 9: Samtools Stats
  SAMTOOLS_STATS(SAMTOOLS_SORT.out.sorted_file)

  // Step 10: Merge BAM files
  // Merge techincal replicates of sample replicates (if tech reps exist). 
  // BAM files for all libraries from same techincal sample replicate. 
  // i.e., merge multiple lanes per sample or resequenceing of 1 sample. 
  // see: https://github.com/nf-core/chipseq/blob/1.2.2/docs/usage.md#multiple-runs-of-the-same-library

  ch_sort_bam_merge = SAMTOOLS_SORT.out.sorted_file
                      .map { it -> [ it[0].split('_')[0..-2].join('_'), it[1] ] }
                      .groupTuple(by: [0])
                      .map { it ->  [ it[0], it[1].flatten() ] }
  // The design script adds 2 fields to sample names: R (replicate), and T (treatment), which are delimited by '_'.
  // The first step, splits off the T identifier, and groups all remaining samples by the new ID. 
  // This allows all samples that are techincal replicataes, i.e., share the same 
  // sampleID and replicate ID, to be joined and then merged in the next step. 

  PICARD_MERGESAMFILES(ch_sort_bam_merge)

  // Step 11: Mark Duplicates
  PICARD_MARKDUPLICATES(PICARD_MERGESAMFILES.out.bam)

  // Step 12: Samtools Stats
  SAMTOOLS_STATS_MD(PICARD_MARKDUPLICATES.out.dedup_bam)

  // Step 13: Samtools Mergebam Filter
  SAMTOOLS_MERGEBAM_FILTER(PICARD_MARKDUPLICATES.out.dedup_bam, MAKE_GENOME_FILTER.out.bed)
  // Note: genome filter file is generic and used by all samples. 

  // JSON files required by BAMTools for alignment filtering
  if (params.read_type == 'SE'){
    ch_bamtools_filter_config = file(params.bamtools_filter_se_config, checkIfExists: true)
  } else {
    ch_bamtools_filter_config = file(params.bamtools_filter_pe_config, checkIfExists: true)
  }

  // Step 14: Bamtools Filter
  BAMTOOLS_FILTER(SAMTOOLS_MERGEBAM_FILTER.out.bam, ch_bamtools_filter_config) 

  // Step 15: Samtools Stats
  SAMTOOLS_STATS_BF(BAMTOOLS_FILTER.out.bam)

  if (params.read_type == 'SE'){

    filtered_sorted_bam = BAMTOOLS_FILTER.out.bam
    // Note: the output BAM from the preceding step was coordinate sorted in step 8, 
    //       and the sort was maintained in the optional merge step and in markduplicates step.

  } else {
    // Step 16: Samtools Name Sort
    NAME_SORT(BAMTOOLS_FILTER.out.bam, '-n -O bam', 'bam')
    // Name sorting is required to remove orphaned singletons in PE data. 

    // Step 17: Remove singleton reads from paired-end BAM file
    BAMPE_RM_ORPHAN(NAME_SORT.out.sorted_file)
    
    // Step 18 : Samtools Pair Sort
    PAIR_SORT(BAMPE_RM_ORPHAN.out.bam, '-O bam', 'bam')
    // Coordinate sorting must be used for next steps. 

    filtered_sorted_bam = PAIR_SORT.out.sorted_file

  }

  // Step 19 : Samtools Stats
  SAMTOOLS_STATS_FILTERED(filtered_sorted_bam)

  // Step 20 : Preseq
  PRESEQ(PICARD_MARKDUPLICATES.out.dedup_bam) 
  //Note: preseq package is aimed at predicting and estimating the complexity of a genomic sequencing library

  // Step 21 : Collect Multiple Metrics

  SAMTOOLS_INDEX(filtered_sorted_bam)

  PICARD_COLLECTMULTIPLEMETRICS(filtered_sorted_bam)

  // Step 22 : Bedtools Genome Coverage
  BEDTOOLS_GENOMECOV(filtered_sorted_bam.join(SAMTOOLS_STATS_FILTERED.out.flagstat))

  // Step 23 : USCS Bedgraph to bigwig
  UCSC_BEDGRAPHTOBIGWIG(BEDTOOLS_GENOMECOV.out.bedgraph, MAKE_GENOME_FILTER.out.sizes) 
  // Note: genome filter is a generic file used for all samples. 

  // Step 24 : Deeptools Compute matrix
  DEEPTOOLS_COMPUTEMATRIX(UCSC_BEDGRAPHTOBIGWIG.out.bigwig, ch_gene_bed)
  // Note: ch_gene_bed is a generic file used for all samples. 

  // Step 25 : Deeptools Plot Profile
  DEEPTOOLS_PLOTPROFILE(DEEPTOOLS_COMPUTEMATRIX.out.matrix)

  // Step 26 : Deeptools Plot Heatmap
  DEEPTOOLS_PLOTHEATMAP(DEEPTOOLS_COMPUTEMATRIX.out.matrix)

  // Step 27 : Phantompeakqualtools
  PHANTOMPEAKQUALTOOLS(filtered_sorted_bam)

  // Step 28 : Multiqc Custom Phantompeakqualtools
  mcp_ch = PHANTOMPEAKQUALTOOLS.out.spp.join(PHANTOMPEAKQUALTOOLS.out.rdata, by: [0])
  MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS(mcp_ch, ch_spp_nsc_header, ch_spp_rsc_header, ch_spp_correlation_header)

  // Create channel linking IP bams with control bams  
  ch_genome_bam_bai = filtered_sorted_bam.join(SAMTOOLS_INDEX.out.bai)
                      .map{it -> [it[0], [it[1], it[2]]]}
                      // next step requires a tuple with [sampleID, [bam, bai]]

  ch_genome_bam_bai = ch_genome_bam_bai
                      .combine(ch_genome_bam_bai)
  // this combine step genenerates pairs of samples, which are then refined in the next step.

  ch_group_bam = control_ch
                    .combine(ch_genome_bam_bai )
                    .filter { it[0] == it[5] && it[1] == it[7] }
                    .join(SAMTOOLS_STATS_FILTERED.out.flagstat)
                    .map { it ->  it[2..-1] }
  // Generate combinations of all study design objects:
  // [SPT5_T0_R1, SPT5_INPUT_R1, SPT5, true, true]
  // with all combined bams from the 'combine' above:
  // [SPT5_T0_R1, [/.../SPT5_T0_R1.mLb.clN.sorted.bam, /.../SPT5_T0_R1.mLb.clN.sorted.bam.bai], SPT5_INPUT_R1, [/.../SPT5_INPUT_R2.mLb.clN.sorted.bam, /.../SPT5_INPUT_R2.mLb.clN.sorted.bam.bai]]
  // the combinations between design and all pairs, has combinations that are not relavent to the study design. Therefore, the combinations are filtered to cases where: 
  // it[0] == it[5] (e.g., SPT5_T0_R1 == SPT5_T0_R1) AND it[1] == it[7] (e.g., SPT5_INPUT_R1 == SPT5_INPUT_R1)
  // it then adjust the output tuple to remove the extra sample IDs: 
  // [SPT5, true, true, SPT5_T0_R2, [/../SPT5_T0_R2.mLb.clN.sorted.bam, /../SPT5_T0_R2.mLb.clN.sorted.bam.bai], SPT5_INPUT_R2, [/../SPT5_INPUT_R2.mLb.clN.sorted.bam, /../SPT5_INPUT_R2.mLb.clN.sorted.bam.bai], /../SPT5_T0_R2.mLb.clN.sorted.bam.flagstat]

  // Step 29 : Deeptools plotFingerprint
  DEEPTOOLS_PLOTFINGERPRINT(ch_group_bam) 

  // Step 30 : Call peaks with MACS2 
  PEAK_CALLING_CHIPSEQ(ch_group_bam, ch_peak_count_header, ch_frip_score_header)
  // Note: ch_peak_count_header is a generic file used for all samples. ch_frip_score_header is a generic file used for all samples. 

  // Step 31 : Calculate FRiP score
  frip_input = ch_group_bam
              .map{it -> [it[3], it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7]]}
              .join(PEAK_CALLING_CHIPSEQ.out.peak)
              .map{it -> it[1..-1]}
  // 'ch_group_bam' is indexed on antibody. peak calling is indexed on the IP sample. 
  // This map adjusts the tuple to put IP in the index position
  // Joins 'ch_group_bam' to the peak file by IP sample ID,
  // and then readjusts the tuple to place antibody in the index position. 

  FRIP_SCORE(frip_input, ch_peak_count_header, ch_frip_score_header)
  // Note: ch_peak_count_header is a generic file used for all samples. ch_frip_score_header is a generic file used for all samples. 

  // Step 32 : Homer Annotate Peaks
  HOMER_ANNOTATEPEAKS(PEAK_CALLING_CHIPSEQ.out.ip_control_peak, ch_fasta, ch_gtf)

  // Step 33 : Plot Macs2 QC
  PLOT_MACS2_QC(PEAK_CALLING_CHIPSEQ.out.peak.collect{ it[-1] })
  // Note: *collect{ it[-1] } collects all peak files, and passes those to the module.  

  // Step 34 : Plot Homer Annotate Peaks
  PLOT_HOMER_ANNOTATEPEAKS(HOMER_ANNOTATEPEAKS.out.txt.collect{ it[-1] }, ch_peak_annotation_header, '_peaks.annotatePeaks.txt')
  // Note: *collect{ it[-1] } collects all peak files, and passes those to the module.  

  // Step 35 : Consensus peaks across samples, create boolean filtering file, SAF file
  
  // Create channel for CONSENSUS PEAKS ANALYSIS 
  // Group by antibody from this point and carry forward boolean variables

  ch_macs_consensus = PEAK_CALLING_CHIPSEQ.out.ip_control_peak
                      .map { it ->  [ it[0], it[1], it[2], it[-1] ] }
                      .groupTuple()
                      .map { it ->  [ it[0], it[1][0], it[2][0], it[3].toSorted( { a, b -> a.getName() <=> b.getName() } ) ] }
  // Note: re-order the output tuple from PEAK_CALLING_CHIPSEQ:
  // [SPT5, true, true, SPT5_T15_R2, SPT5_INPUT_R2, /.../SPT5_T15_R2_peaks.broadPeak]
  // to remove the case and control sample IDs:
  // [SPT5, true, true, /.../SPT5_T15_R1_peaks.broadPeak]
  // Then group by antibody. Map: keep only the first index position of replicatesExist, multipleGroups
  // as the remaining array for those are duplicate values. sort the broadpeak file array by file name. 

  MACS2_CONSENSUS(ch_macs_consensus)
  // Note: this step will not run when replicatesExist || multipleGroups are false. 
  //       Subequently all steps beyond this point will not run as they rely on output from this step. 

  // Step 36 : Consensus peaks annotation
  CONSENSUS_PEAKS_ANNOTATE(MACS2_CONSENSUS.out.bed, ch_fasta, ch_gtf)
  // Note: ch_fasta and ch_gtf are generic files and shared by all samples. 

  // Step 37 : Annotate boolean peaks
  ANNOTATE_BOOLEAN_PEAKS(MACS2_CONSENSUS.out.boolean_txt.join(CONSENSUS_PEAKS_ANNOTATE.out.txt))
  
  // Get BAM and SAF files for each antibody

  ch_group_bam                                                 // [antibody, replicatesExist, multipleGroups, sample_id, [bam, bai], control_id, [bam, bai], sample_id bam.flagstat] 
      .map { it -> [ it[3], [ it[0], it[1], it[2] ] ] }        // [sample_id, [antibody, replicatesExist, multipleGroups]] 
      .join(filtered_sorted_bam)                               // [sample_id, [antibody, replicatesExist, multipleGroups], final filtered sample_id indexed bam]
      .map { it -> [ it[1][0], it[1][1], it[1][2], it[2] ] }   // [antibody, replicatesExist, multipleGroups, OR sample_id bam]
      .groupTuple()
      .map { it -> [ it[0], it[1][0], it[2][0], it[3].flatten().sort() ] } // [antibody, replicatesExist, multipleGroups, [OR sample_id1 R1 bam, OR sample_id1 R2 bam, OR sample_id2 R1 bam, OR sample_id2 R2 bam]]
      .join(MACS2_CONSENSUS.out.saf)  // [antibody, replicatesExist, multipleGroups, [OR sample_id1 R1 bam, OR sample_id1 R2 bam, OR sample_id2 R1 bam, OR sample_id2 R2 bam], SAF]
      .set { ch_group_bam }

  // Step 38 : Count reads in consensus peaks with featureCounts
  SUBREAD_FEATURECOUNTS(ch_group_bam)

  // Step 39 : Differential analysis with DESeq2
  DESEQ2_QC(SUBREAD_FEATURECOUNTS.out.counts, ch_deseq2_pca_header, ch_deseq2_clustering_header)
  // note: ch_deseq2_pca_header, ch_deseq2_clustering_header are generic files used for all samples. 

  // Create channels for multi input files
  ch_multiqc_files = Channel.empty()
  
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(TRIM_GALORE.out.trim_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(TRIM_GALORE.out.trimmed_fastqc.collect{it[1]}.ifEmpty([]))

  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.flagstat.collect{it[1]}.ifEmpty([]))  
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.idxstat.collect{it[1]}.ifEmpty([]))  
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{it[1]}.ifEmpty([]))  
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_MD.out.flagstat.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_MD.out.idxstat.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_MD.out.stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_FILTERED.out.flagstat.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_FILTERED.out.idxstat.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_FILTERED.out.stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.dedup_metrics.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics.collect{it[1]}.ifEmpty([]))

  ch_multiqc_files = ch_multiqc_files.mix(FRIP_SCORE.out.tsv.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PLOT_HOMER_ANNOTATEPEAKS.out.tsv.collect())
  ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.pca_multiqc.collect())
  ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.dists_multiqc.collect())

  ch_multiqc_files = ch_multiqc_files.mix(PRESEQ.out.txt.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTFINGERPRINT.out.raw.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTPROFILE.out.table.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PHANTOMPEAKQUALTOOLS.out.spp.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.nsc.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.rsc.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.correlation.collect{it[1]}.ifEmpty([]))
   

  // Step 41 : MultiQC  
  MULTIQC (
      ch_multiqc_files.collect()
  )

}
