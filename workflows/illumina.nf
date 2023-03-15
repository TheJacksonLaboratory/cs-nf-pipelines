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
include {PICARD_SORTSAM as LUMPY_SORTSAM;
         PICARD_SORTSAM as LUMPY_SORTSAM_DISCORDANT;
         PICARD_SORTSAM as LUMPY_SORTSAM_SPLIT} from "${projectDir}/modules/picard/picard_sortsam"
include {LUMPY_EXTRACT_SPLITS} from "${projectDir}/modules/lumpy/lumpy_extract_splits"
include {LUMPY_CALL_SV} from "${projectDir}/modules/lumpy/lumpy_call_sv"
include {REHEADER_VCF as REHEADER_LUMPY;
         REHEADER_VCF as REHEADER_BREAKDANCER} from "${projectDir}/modules/utility_modules/reheader_vcf"
include {BREAKDANCER_CALL} from "${projectDir}/modules/breakdancer/breakdancer_call"
include {BREAKDANCER_SV_TO_VCF} from "${projectDir}/modules/breakdancer/breakdancer_sv_to_vcf"

workflow ILLUMINA {
    params.fasta = params.genome ? params.genomes[params.genome].fasta ?: null : null
    ch_fasta = Channel.fromPath(params.fasta)
    ch_bwa_index = params.bwa_index ? Channel.fromPath(params.bwa_index) : null
    ch_fastq1 = params.fastq1 ? Channel.fromPath(params.fastq1) : null
    ch_fastq2 = params.fastq2 ? Channel.fromPath(params.fastq2) : null
    ch_sampleID = params.names ? Channel.value(params.names) : null
    ch_bam = params.bam ? Channel.fromPath(params.bam) : null
    ch_bwa_index = params.bwa_index ? Channel.fromPath(params.bwa_index) : null

    // Prepare reads channel

    if (!params.fastq2 && !params.bam) {
        fq_reads = ch_sampleID.concat(ch_fastq1)
                            .collect()
                            .map { it -> tuple(it[0], it[1])}
    }

    else if (params.fastq1 && params.fastq2) {
        fq_reads = ch_sampleID.concat(ch_fastq1, ch_fastq2)
                            .collect()
                            .map { it -> tuple(it[0], tuple(it[1], it[2]))}
    }

    else {
        fq_reads = null
        pre_bam = ch_sampleID.concat(ch_bam)
                             .collect()
                             .map { it -> tuple(it[0], it[1])}
    }

    // Step 0: Generate reference index if neccesary
    if(!ch_bwa_index) {
        BWA_INDEX(ch_fasta)
        ch_bwa_index = BWA_INDEX.out.bwa_index
    }

    // Index reference fasta
    SAMTOOLS_FAIDX(ch_fasta)
    
    // ** Optional mapping steps when input are FASTQ files
    if (params.fastq1) {
        // Get read groups ID from FASTQ file
        READ_GROUPS(fq_reads)

        // Map reads to reference
        bwa_mem_input = fq_reads.join(READ_GROUPS.out.read_groups)
        BWA_MEM(bwa_mem_input, SAMTOOLS_FAIDX.out.fasta_fai, ch_bwa_index)

        // Sort and compress to BAM
        SAMTOOLS_SORT(BWA_MEM.out.sam)

        ch_bam_undup = SAMTOOLS_SORT.out.bam
    }
    else {
        ch_bam_undup = ch_bam
    }

    // Remove optical duplicates from alignment
    GATK_MARK_DUPLICATES(ch_bam_undup)

    // Quantify insert sizes
    SAMTOOLS_STATS(GATK_MARK_DUPLICATES.out.bam_and_index)

    // Prep BAM for Lumpy (Map clipped reads, read group info, extract discordant alignments)
    LUMPY_PREP(GATK_MARK_DUPLICATES.out.bam_and_index)

    // Sort prepped LUMPY bams
    LUMPY_SORTSAM(LUMPY_PREP.out.bam_bwa_lumpy, "alignBWA_lumpySort")
    LUMPY_SORTSAM_DISCORDANT(LUMPY_PREP.out.dis_unsorted_bam, "alignBWA_lumpySort_discordant")

    // Extract split reads
    LUMPY_EXTRACT_SPLITS(LUMPY_PREP.out.bam_bwa_lumpy)
    LUMPY_SORTSAM_SPLIT(LUMPY_EXTRACT_SPLITS.out.bam_bwa_lumpy, "alignBWA_lumpySort_splits")

    // Call SV with Lumpy
    lumpy_input = LUMPY_SORTSAM.out.sorted_bam.join(LUMPY_SORTSAM_SPLIT.out.sorted_bam).join(LUMPY_SORTSAM_DISCORDANT.out.sorted_bam)
    LUMPY_CALL_SV(lumpy_input)
    REHEADER_LUMPY(LUMPY_CALL_SV.out.lumpy_vcf, "lumpy")

    // * Breakdancer

    // Call SV with Breakdancer
    BREAKDANCER_CALL(GATK_MARK_DUPLICATES.out.bam_and_index)

    // Convert Breakdancer SV to VCF
    breakdancer_vcf_input = GATK_MARK_DUPLICATES.out.bam_and_index.join(BREAKDANCER_CALL.out.breakdancer_sv)
    BREAKDANCER_SV_TO_VCF(breakdancer_vcf_input)
    REHEADER_BREAKDANCER(BREAKDANCER_SV_TO_VCF.out.breakdancer_vcf, "breakdancer")
}