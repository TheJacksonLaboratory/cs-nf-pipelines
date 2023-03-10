#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//include {param_log} from 
include {BWA_INDEX} from "${projectDir}/modules/bwa/bwa_index"
include {SAMTOOLS_FAIDX} from "${projectDir}/modules/samtools/samtools_faidx"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {SAMTOOLS_SORT} from "${projectDir}/modules/samtools/samtools_sort"
include {GATK_MARK_DUPLICATES} from "${projectDir}/modules/gatk/gatk_mark_duplicates"
include {SAMTOOLS_STATS} from "${projectDir}/modules/samtools/samtools_stats"
include {LUMPY_PREP} from "${projectDir}/modules/lumpy/lumpy_prep"
include {PICARD_SORTSAM as LUMPY_SORTSAM,
         PICARD_SORTSAM as LUMPY_SORTSAM_DISCORDANT,
         PICARD_SORTSAM as LUMPY_SORTSAM_SPLIT} from "${projectDir}/modules/picard/picard_sortsam"
include {LUMPY_EXTRACT_SPLITS} from "${projectDir}/modules/lumpy/lumpy_extract_splits"

workflow ILLUMINA {
    params.fasta = params.genome ? params.genomes[params.genome].fasta ?: null : null
    ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : "null"
    ch_fastq1 = params.fastq1 ? Channel.value(file(params.fastq1)) : null
    ch_fastq2 = params.fastq2 ? Channel.value(file(params.fastq2)) : null
    
    // Prepare reads channel

    if (!params.fastq2 && !params.bam) {
        fq_reads = tuple(val(params.names), file(params.fastq1))
    }

    else if (!params.bam) {
        fq_reads = tuple(val(params.names), file(params.fastq1), file(params.fastq2))
    }

    else {
        fq_reads = null
        pre_bam = tuple(val(params.names), file(params.bam))
    }

    // Step 0: Generate reference index if neccesary
    if(!params.bwa_index || params.genome) {
        BWA_INDEX(params.fasta)
        params.bwa_index = BWA_INDEX.out.bwa_index
    }

    // Index reference fasta
    SAMTOOLS_FAIDX(params.fasta)

    // ** Optional mapping steps when input are FASTQ files
    if (params.fastq1) {
        // Get read groups ID from FASTQ file
        READ_GROUPS(fq_reads)

        // Map reads to reference
        bwa_mem_input = fq_reads.join(READ_GROUPS.out.read_groups)
        BWA_MEM(bwa_mem_input, params.bwa_index)

        // Sort and compress to BAM
        SAMTOOLS_SORT(BWA_MEM.out.sam)

        ch_bam_undup = SAMTOOLS_SORT.out.bam
    }
    else {
        ch_bam_undup = params.bam ? Channel.value(file(params.bam)) : null
    }

    // Remove optical duplicates from alignment
    GATK_MARK_DUPLICATES(ch_bam_undup)

    // Quantify insert sizes
    SAMTOOLS_STATS(GATK_MARK_DUPLICATES.out.bam_and_index)

    // Prep BAM for Lumpy (Map clipped reads, read group info, extract discordant alignments)
    LUMPY_PREP(GATK_MARK_DUPLICATES.out.bam_and_index)

    // Sort prepped LUMPY bams
    LUMPY_SORTSAM(LUMPY_PREP.out.bam_bwa_lumpy)
    LUMPY_SORTSAM_DISCORDANT(LUMPY_PREP.out.dis_unsorted_bam)

    // Extract split reads
    LUMPY_EXTRACT_SPLITS(LUMPY_PREP.out.bam_bwa_lumpy)
}