#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {QUALITY_STATISTICS} from "${projectDir}/modules/utility_modules/quality_stats"
include {READ_GROUPS as READ_GROUPS_HUMAN;
         READ_GROUPS as READ_GROUPS_MOUSE} from "${projectDir}/modules/utility_modules/read_groups"
include {FASTQ_PAIR} from "${projectDir}/modules/fastq-tools/fastq-pair"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {XENOME_CLASSIFY} from "${projectDir}/modules/xenome/xenome"
include {FASTQ_SORT as FASTQ_SORT_HUMAN;
         FASTQ_SORT as FASTQ_SORT_MOUSE} from "${projectDir}/modules/fastq-tools/fastq-sort"
include {RSEM_ALIGNMENT_EXPRESSION as RSEM_ALIGNMENT_EXPRESSION_HUMAN;
         RSEM_ALIGNMENT_EXPRESSION as RSEM_ALIGNMENT_EXPRESSION_MOUSE} from "${projectDir}/modules/rsem/rsem_alignment_expression"
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
    // Step 1: Qual_Stat, Get read group information, Run Xenome
    QUALITY_STATISTICS(read_ch)

    FASTQ_PAIR(QUALITY_STATISTICS.out.trimmed_fastq)

    // QC is assess on all reads. Mouse/human is irrelevant here. 
    FASTQC(QUALITY_STATISTICS.out.trimmed_fastq)

    // Xenome Classification
    XENOME_CLASSIFY(FASTQ_PAIR.out.paired_fastq)

    // Xenome Read Sort
    FASTQ_SORT_HUMAN(XENOME_CLASSIFY.out.xenome_fastq, 'human')
    FASTQ_SORT_MOUSE(XENOME_CLASSIFY.out.xenome_mouse_fastq, 'mouse')    

    human_read = FASTQ_SORT_HUMAN.out.sorted_fastq
    .map{it -> tuple(it[0]+'_human', it[1])}

    mouse_reads = FASTQ_SORT_MOUSE.out.sorted_fastq
    .map{it -> tuple(it[0]+'_mouse', it[1])}

    // Step 2: RSEM Human and Stats: 

    if (params.rsem_aligner == "bowtie2") {
    rsem_ref_files = file("${params.rsem_ref_files_human}/bowtie2/*")
    }
    else if (params.rsem_aligner == "star") {
    rsem_ref_files = file("${params.rsem_ref_files_human}/STAR/${params.rsem_star_prefix_human}/*")
    }
    else error "${params.rsem_aligner_human} is not valid, use 'bowtie2' or 'star'"

    RSEM_ALIGNMENT_EXPRESSION_HUMAN(human_read, rsem_ref_files, params.rsem_ref_prefix_human)
    
    // Picard Alignment Metrics
    READ_GROUPS_HUMAN(human_read, "picard")

    add_replace_groups_human = READ_GROUPS_HUMAN.out.read_groups.join(RSEM_ALIGNMENT_EXPRESSION_HUMAN.out.bam)
    PICARD_ADDORREPLACEREADGROUPS_HUMAN(add_replace_groups_human)

    PICARD_REORDERSAM_HUMAN(PICARD_ADDORREPLACEREADGROUPS_HUMAN.out.bam, params.picard_dict_human)

    // Picard Alignment Metrics
    PICARD_SORTSAM_HUMAN(PICARD_REORDERSAM_HUMAN.out.bam)
    PICARD_COLLECTRNASEQMETRICS_HUMAN(PICARD_SORTSAM_HUMAN.out.bam, params.ref_flat_human, params.ribo_intervals_human)

    // Step 3 RSEM Mouse and Stats:
    if (params.rsem_aligner == "bowtie2") {
    rsem_ref_files = file("${params.rsem_ref_files_mouse}/bowtie2/*")
    }
    else if (params.rsem_aligner == "star") {
    rsem_ref_files = file("${params.rsem_ref_files_mouse}/STAR/${params.rsem_star_prefix_mouse}/*")
    }
    else error "${params.rsem_aligner_mouse} is not valid, use 'bowtie2' or 'star'"

    RSEM_ALIGNMENT_EXPRESSION_MOUSE(mouse_reads, rsem_ref_files, params.rsem_ref_prefix_mouse)
    
    // Step 4: Picard Alignment Metrics
    READ_GROUPS_MOUSE(mouse_reads, "picard")

    add_replace_groups_mouse = READ_GROUPS_MOUSE.out.read_groups.join(RSEM_ALIGNMENT_EXPRESSION_MOUSE.out.bam)
    PICARD_ADDORREPLACEREADGROUPS_MOUSE(add_replace_groups_mouse)

    PICARD_REORDERSAM_MOUSE(PICARD_ADDORREPLACEREADGROUPS_MOUSE.out.bam, params.picard_dict_mouse)

    // Step 5: Picard Alignment Metrics
    PICARD_SORTSAM_MOUSE(PICARD_REORDERSAM_MOUSE.out.bam)
    PICARD_COLLECTRNASEQMETRICS_MOUSE(PICARD_SORTSAM_MOUSE.out.bam, params.ref_flat_mouse, params.ribo_intervals_mouse)


    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(QUALITY_STATISTICS.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(XENOME_CLASSIFY.out.xenome_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(RSEM_ALIGNMENT_EXPRESSION_HUMAN.out.rsem_cnt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTRNASEQMETRICS_HUMAN.out.picard_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(RSEM_ALIGNMENT_EXPRESSION_MOUSE.out.rsem_cnt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTRNASEQMETRICS_MOUSE.out.picard_metrics.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )

}