#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {FASTP} from "${projectDir}/modules/fastp/fastp"
include {READ_GROUPS as READ_GROUPS_HUMAN;
         READ_GROUPS as READ_GROUPS_MOUSE} from "${projectDir}/modules/utility_modules/read_groups"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {GET_READ_LENGTH} from "${projectDir}/modules/utility_modules/get_read_length"
include {CHECK_STRANDEDNESS} from "${projectDir}/modules/python/python_check_strandedness"
include {XENGSORT_INDEX} from "${projectDir}/modules/xengsort/xengsort_index"
include {XENGSORT_CLASSIFY} from "${projectDir}/modules/xengsort/xengsort_classify"
// include {GZIP as GZIP_HUMAN;
//          GZIP as GZIP_MOUSE} from "${projectDir}/modules/utility_modules/gzip"
include {RSEM_ALIGNMENT_EXPRESSION as RSEM_ALIGNMENT_EXPRESSION_HUMAN;
         RSEM_ALIGNMENT_EXPRESSION as RSEM_ALIGNMENT_EXPRESSION_MOUSE} from "${projectDir}/modules/rsem/rsem_alignment_expression"
include {LYMPHOMA_CLASSIFIER} from "${projectDir}/modules/python/python_lymphoma_classifier"
include {PICARD_ADDORREPLACEREADGROUPS as PICARD_ADDORREPLACEREADGROUPS_HUMAN;
         PICARD_ADDORREPLACEREADGROUPS as PICARD_ADDORREPLACEREADGROUPS_MOUSE} from "${projectDir}/modules/picard/picard_addorreplacereadgroups"
include {PICARD_REORDERSAM as PICARD_REORDERSAM_HUMAN;
         PICARD_REORDERSAM as PICARD_REORDERSAM_MOUSE} from "${projectDir}/modules/picard/picard_reordersam"
include {PICARD_SORTSAM as PICARD_SORTSAM_HUMAN;
         PICARD_SORTSAM as PICARD_SORTSAM_MOUSE} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_COLLECTRNASEQMETRICS as PICARD_COLLECTRNASEQMETRICS_HUMAN;
         PICARD_COLLECTRNASEQMETRICS as PICARD_COLLECTRNASEQMETRICS_MOUSE} from "${projectDir}/modules/picard/picard_collectrnaseqmetrics"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

workflow PDX_RNASEQ {

    take:
        read_ch

    main:
    // Step 1: Read trim, Get read group information, Run xengsort
    FASTP(read_ch)
    
    GET_READ_LENGTH(read_ch)
    
    if (params.read_type == 'PE') {
      xengsort_input = FASTP.out.trimmed_fastq
    } else {
      xengsort_input = FASTP.out.trimmed_fastq
    }

    // QC is assess on all reads. Mouse/human is irrelevant here. 
    FASTQC(FASTP.out.trimmed_fastq)

    CHECK_STRANDEDNESS(FASTP.out.trimmed_fastq)

    // Generate Xengsort Index if needed
    if (params.xengsort_idx_path) {
        xengsort_index = params.xengsort_idx_path
    } else {
        XENGSORT_INDEX(params.xengsort_host_fasta, params.ref_fa)
        xengsort_index = XENGSORT_INDEX.out.xengsort_index
    }

    // Xengsort Classification
    XENGSORT_CLASSIFY(xengsort_index, xengsort_input) 

    human_reads = XENGSORT_CLASSIFY.out.xengsort_human_fastq
                  .join(CHECK_STRANDEDNESS.out.strand_setting)
                  .join(GET_READ_LENGTH.out.read_length)
                  .map{it -> tuple(it[0]+'_human', it[1], it[2], it[3])}

    mouse_reads = XENGSORT_CLASSIFY.out.xengsort_mouse_fastq
                  .join(CHECK_STRANDEDNESS.out.strand_setting)
                  .join(GET_READ_LENGTH.out.read_length)
                  .map{it -> tuple(it[0]+'_mouse', it[1], it[2], it[3])}

    // Step 2: RSEM Human and Stats: 

    RSEM_ALIGNMENT_EXPRESSION_HUMAN(human_reads, params.rsem_ref_files_human, params.rsem_star_prefix_human, params.rsem_ref_prefix_human)
    
    LYMPHOMA_CLASSIFIER(RSEM_ALIGNMENT_EXPRESSION_HUMAN.out.rsem_genes)

    // Picard Alignment Metrics
    READ_GROUPS_HUMAN(human_reads.map{it -> tuple(it[0], it[1])}, "picard")

    add_replace_groups_human = READ_GROUPS_HUMAN.out.read_groups.join(RSEM_ALIGNMENT_EXPRESSION_HUMAN.out.bam)
    PICARD_ADDORREPLACEREADGROUPS_HUMAN(add_replace_groups_human)

    PICARD_REORDERSAM_HUMAN(PICARD_ADDORREPLACEREADGROUPS_HUMAN.out.bam, params.picard_dict_human)

    // Picard Alignment Metrics
    PICARD_SORTSAM_HUMAN(PICARD_REORDERSAM_HUMAN.out.bam, 'coordinate')

    human_qc_input = PICARD_SORTSAM_HUMAN.out.bam.join(human_reads)
                     .map{it -> [it[0], it[1], it[3]]}
                     
    PICARD_COLLECTRNASEQMETRICS_HUMAN(human_qc_input, params.ref_flat_human, params.ribo_intervals_human)

    // Step 3 RSEM Mouse and Stats:

    RSEM_ALIGNMENT_EXPRESSION_MOUSE(mouse_reads, params.rsem_ref_files_mouse, params.rsem_star_prefix_mouse, params.rsem_ref_prefix_mouse)
    
    // Step 4: Picard Alignment Metrics
    READ_GROUPS_MOUSE(mouse_reads.map{it -> tuple(it[0], it[1])}, "picard")

    add_replace_groups_mouse = READ_GROUPS_MOUSE.out.read_groups.join(RSEM_ALIGNMENT_EXPRESSION_MOUSE.out.bam)
    PICARD_ADDORREPLACEREADGROUPS_MOUSE(add_replace_groups_mouse)

    PICARD_REORDERSAM_MOUSE(PICARD_ADDORREPLACEREADGROUPS_MOUSE.out.bam, params.picard_dict_mouse)

    // Step 5: Picard Alignment Metrics
    PICARD_SORTSAM_MOUSE(PICARD_REORDERSAM_MOUSE.out.bam, 'coordinate')

    mouse_qc_input = PICARD_SORTSAM_MOUSE.out.bam.join(mouse_reads)
                     .map{it -> [it[0], it[1], it[3]]}
 
    PICARD_COLLECTRNASEQMETRICS_MOUSE(mouse_qc_input, params.ref_flat_mouse, params.ribo_intervals_mouse)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.quality_json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(XENGSORT_CLASSIFY.out.xengsort_log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(RSEM_ALIGNMENT_EXPRESSION_HUMAN.out.rsem_cnt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(RSEM_ALIGNMENT_EXPRESSION_HUMAN.out.star_log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTRNASEQMETRICS_HUMAN.out.picard_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(RSEM_ALIGNMENT_EXPRESSION_MOUSE.out.rsem_cnt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(RSEM_ALIGNMENT_EXPRESSION_MOUSE.out.star_log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTRNASEQMETRICS_MOUSE.out.picard_metrics.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )

} 
