#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
//include {help} from '../bin/help/rnaseq'

include {NANOSTAT as NANOSTAT_PREFILT;
         NANOSTAT as NANOSTAT_POSTFILT} from "${projectDir}/modules/nanostat/nanostat"
include {PORECHOP} from "${projectDir}/modules/porechop/porechop"
include {NANOQC} from "${projectDir}/modules/nanoqc/nanoqc"
include {NANOFILT} from "${projectDir}/modules/nanofilt/nanofilt"
include {MINIMAP2_INDEX} from "${projectDir}/modules/minimap/minimap2_index"
include {MINIMAP2_MAP_ONT} from "${projectDir}/modules/minimap/minimap2_map_ont"
include {SAMTOOLS_SORT} from "${projectDir}/modules/samtools/samtools_sort"
include {SURVIVOR_MERGE} from "${projectDir}/modules/survivor/survivor_merge"
include {SURVIVOR_VCF_TO_TABLE} from "${projectDir}/modules/survivor/survivor_vcf_to_table"
include {SURVIVOR_SUMMARY} from "${projectDir}/modules/survivor/survivor_summary"
include {SURVIVOR_TO_BED} from "${projectDir}/modules/survivor/survivor_to_bed"
include {SURVIVOR_BED_INTERSECT} from "${projectDir}/modules/survivor/survivor_bed_intersect"
include {SURVIVOR_ANNOTATION} from "${projectDir}/modules/survivor/survivor_annotation"
include {SURVIVOR_INEXON} from "${projectDir}/modules/survivor/survivor_inexon"

// log paramater info
//param_log()

workflow ONT {
    params.fasta = params.genome ? params.genomes[params.genome].fasta ?: null : null
    ch_fasta = Channel.fromPath(params.fasta)
    ch_fastq1 = params.fastq1 ? Channel.fromPath(params.fastq1) : null
    ch_sampleID = params.names ? Channel.value(params.names) : null
    ch_bam = params.bam ? Channel.fromPath(params.bam) : null

    if (params.fastq1 && !params.bam) {
        fq_reads = ch_sampleID.concat(ch_fastq1)
                            .collect()
                            .map { it -> tuple(it[0], it[1])}
    }

    else {
        fq_reads = null
        pre_bam = ch_sampleID.concat(ch_bam)
                             .collect()
                             .map { it -> tuple(it[0], it[1])}
    }

    // ** Optional mapping steps when input is a FASTQ file
    if (params.fastq1) {
        
        NANOSTAT_PREFILT(fq_reads)

        PORECHOP(fq_reads)

        NANOQC(PORECHOP.out.porechop_fastq)

        NANOFILT(PORECHOP.out.porechop_fastq)

        NANOSTAT_POSTFILT(NANOFILT.out.porechop_nanofilt_fastq)

        // Prepare index
        MINIMAP2_INDEX(ch_fasta)

        // Map reads to indexed genome
        MINIMAP2_MAP_ONT(NANOFILT.out.porechop_nanofilt_fastq, MINIMAP2_INDEX.out.minimap2_index)

        SAMTOOLS_SORT(MINIMAP2_MAP_ONT.out.minimap_sam)
        
        ch_mm2_bam =  SAMTOOLS_SORT.out.bam
    }

    else {
        ch_mm2_bam = pre_bam
    }

}