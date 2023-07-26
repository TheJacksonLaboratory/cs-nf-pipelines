#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/ont"
include {PARAM_LOG} from "${projectDir}/bin/log/ont"
include {NANOSTAT as NANOSTAT_PREFILT;
         NANOSTAT as NANOSTAT_POSTFILT} from "${projectDir}/modules/nanostat/nanostat"
include {PORECHOP} from "${projectDir}/modules/porechop/porechop"
include {NANOQC} from "${projectDir}/modules/nanoqc/nanoqc"
include {NANOFILT} from "${projectDir}/modules/nanofilt/nanofilt"
include {MINIMAP2_INDEX} from "${projectDir}/modules/minimap/minimap2_index"
include {MINIMAP2_MAP_ONT} from "${projectDir}/modules/minimap/minimap2_map_ont"
include {SAMTOOLS_SORT} from "${projectDir}/modules/samtools/samtools_sort"
include {SAMTOOLS_FILTER} from "${projectDir}/modules/samtools/samtools_filter"
include {SNIFFLES} from "${projectDir}/modules/sniffles/sniffles"
include {NANOSV} from "${projectDir}/modules/nanosv/nanosv"
include {PYTHON_PARSE_DEPTHS as PARSE_NANOSV_DEPTHS;
         PYTHON_PARSE_DEPTHS as PARSE_SNIFFLES_DEPTHS} from "${projectDir}/modules/python/python_parse_depths"
include {SURVIVOR_MERGE} from "${projectDir}/modules/survivor/survivor_merge"
include {PYTHON_PARSE_SURVIVOR_IDS} from "${projectDir}/modules/python/python_parse_survivor_ids"
include {VCFTOOLS_FILTER} from "${projectDir}/modules/vcftools/vcftools_filter"
include {SURVIVOR_VCF_TO_TABLE} from "${projectDir}/modules/survivor/survivor_vcf_to_table"
include {SURVIVOR_SUMMARY} from "${projectDir}/modules/survivor/survivor_summary"
include {SURVIVOR_TO_BED} from "${projectDir}/modules/survivor/survivor_to_bed"
include {SURVIVOR_BED_INTERSECT} from "${projectDir}/modules/survivor/survivor_bed_intersect"
include {SURVIVOR_ANNOTATION} from "${projectDir}/modules/survivor/survivor_annotation"
include {R_MERGE_DEPTHS} from "${projectDir}/modules/r/r_merge_depths"
include {SURVIVOR_INEXON} from "${projectDir}/modules/survivor/survivor_inexon"
include {PYTHON_ANNOT_DEPTHS} from "${projectDir}/modules/python/python_annot_depths"
include {PYTHON_ANNOT_ON_TARGET} from "${projectDir}/modules/python/python_annot_on_target"

// log paramater info
//param_log()

PARAM_LOG()

workflow ONT {
    if (params.help){
       help()
       exit 0
    }

    ch_fasta = params.fasta ? Channel.fromPath(params.fasta): null
    ch_fastq1 = params.fastq1 ? Channel.fromPath(params.fastq1) : null
    ch_sampleID = params.sampleID ? Channel.value(params.sampleID) : null
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
        

        // Generate reference index if neccesary
        if(!params.minimap2_index) {
            MINIMAP2_INDEX(ch_fasta)
            ch_minimap2_index = MINIMAP2_INDEX.out.minimap2_index
        }
        else {
            ch_minimap2_index = file("${params.minimap2_index}")
        }        

        // Map reads to indexed genome
        MINIMAP2_MAP_ONT(NANOFILT.out.porechop_nanofilt_fastq, ch_minimap2_index)

        SAMTOOLS_SORT(MINIMAP2_MAP_ONT.out.minimap_sam)
        
        ch_mm2_bam =  SAMTOOLS_SORT.out.bam
    }

    else {
        ch_mm2_bam = pre_bam
    }

    SAMTOOLS_FILTER(ch_mm2_bam)

    SNIFFLES(SAMTOOLS_FILTER.out.bam_and_index)
    PARSE_SNIFFLES_DEPTHS(SNIFFLES.out.sniffles_vcf, "sniffles")

    NANOSV(SAMTOOLS_FILTER.out.bam_and_index)
    PARSE_NANOSV_DEPTHS(NANOSV.out.nanosv_vcf, "nanosv")

    // Join VCFs together by sampleID and run SURVIVOR merge

    survivor_input = NANOSV.out.nanosv_vcf.join(SNIFFLES.out.sniffles_vcf)
                     .map { it -> tuple(it[0], tuple(it[1], it[2]))}
    SURVIVOR_MERGE(survivor_input)
    PYTHON_PARSE_SURVIVOR_IDS(SURVIVOR_MERGE.out.vcf)
    SURVIVOR_VCF_TO_TABLE(SURVIVOR_MERGE.out.vcf)
    SURVIVOR_SUMMARY(SURVIVOR_MERGE.out.vcf)

    bed_prep_input = SURVIVOR_VCF_TO_TABLE.out.annotation.join(SURVIVOR_SUMMARY.out.csv)
    SURVIVOR_TO_BED(bed_prep_input)
    SURVIVOR_BED_INTERSECT(SURVIVOR_TO_BED.out.sv_beds)
    surv_annot_input = SURVIVOR_TO_BED.out.sv_beds.join(SURVIVOR_BED_INTERSECT.out.intersected_beds).join(SURVIVOR_SUMMARY.out.csv).join(SURVIVOR_VCF_TO_TABLE.out.annotation)
    SURVIVOR_ANNOTATION(surv_annot_input)
    
    surv_depths_input = PARSE_NANOSV_DEPTHS.out.csv.join(PARSE_SNIFFLES_DEPTHS.out.csv).join(PYTHON_PARSE_SURVIVOR_IDS.out.csv).join(SURVIVOR_ANNOTATION.out.csv)
    R_MERGE_DEPTHS(surv_depths_input)

    VCFTOOLS_FILTER(SURVIVOR_MERGE.out.vcf)
    surv_inexon_input =  VCFTOOLS_FILTER.out.vcf.join(SURVIVOR_BED_INTERSECT.out.intersected_exons)
    SURVIVOR_INEXON(surv_inexon_input)

    depths_input = SURVIVOR_INEXON.out.vcf.join(R_MERGE_DEPTHS.out.bed)
    PYTHON_ANNOT_DEPTHS(depths_input)

    PYTHON_ANNOT_ON_TARGET(PYTHON_ANNOT_DEPTHS.out.vcf)
}