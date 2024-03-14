#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/pacbio"
include {PARAM_LOG} from "${projectDir}/bin/log/pacbio"
include {PBMM2_INDEX} from "${projectDir}/modules/pbmm2/pbmm2_index"
include {PBMM2_CALL} from "${projectDir}/modules/pbmm2/pbmm2_call"
include {PBSV_DISCOVER} from "${projectDir}/modules/pbsv/pbsv_discover"
include {PBSV_CALL} from "${projectDir}/modules/pbsv/pbsv_call"
include {SNIFFLES} from "${projectDir}/modules/sniffles/sniffles"
include {SURVIVOR_MERGE} from "${projectDir}/modules/survivor/survivor_merge"
include {SURVIVOR_VCF_TO_TABLE} from "${projectDir}/modules/survivor/survivor_vcf_to_table"
include {SURVIVOR_SUMMARY} from "${projectDir}/modules/survivor/survivor_summary"
include {SURVIVOR_TO_BED} from "${projectDir}/modules/survivor/survivor_to_bed"
include {SURVIVOR_BED_INTERSECT} from "${projectDir}/modules/survivor/survivor_bed_intersect"
include {SURVIVOR_ANNOTATION} from "${projectDir}/modules/survivor/survivor_annotation"
include {SURVIVOR_INEXON} from "${projectDir}/modules/survivor/survivor_inexon"

// log parameter info
PARAM_LOG()

workflow PACBIO {

    if (params.help){
       help()
        exit 0
    }

    ch_fasta = params.ref_fa ? Channel.fromPath(params.ref_fa): null
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
        ch_minimap2_index = file("${params.minimap2_index}")

        // Map reads to indexed genome
        PBMM2_CALL(fq_reads, ch_minimap2_index)
        // Map reads to reference

        ch_pbmm2_bam = PBMM2_CALL.out.pbmm2_bam
    }

    else {
        ch_pbmm2_bam = pre_bam
    }

    // Call SV with PBSV
    PBSV_DISCOVER(ch_pbmm2_bam)
    PBSV_CALL(PBSV_DISCOVER.out.pbsv_svsig, ch_fasta)

    // Call SV with sniffles
    SNIFFLES(ch_pbmm2_bam)

    // * Merge callers and annotate results

    // Join VCFs together by sampleID and run SURVIVOR merge

    survivor_input = PBSV_CALL.out.pbsv_vcf.join(SNIFFLES.out.sniffles_vcf)
                     .map { it -> tuple(it[0], tuple(it[1], it[2]))}
    SURVIVOR_MERGE(survivor_input)
    SURVIVOR_VCF_TO_TABLE(SURVIVOR_MERGE.out.vcf)
    SURVIVOR_SUMMARY(SURVIVOR_MERGE.out.vcf)

    bed_prep_input = SURVIVOR_VCF_TO_TABLE.out.annotation.join(SURVIVOR_SUMMARY.out.csv)
    SURVIVOR_TO_BED(bed_prep_input)
    SURVIVOR_BED_INTERSECT(SURVIVOR_TO_BED.out.sv_beds)
    surv_annot_input = SURVIVOR_TO_BED.out.sv_beds.join(SURVIVOR_BED_INTERSECT.out.intersected_beds).join(SURVIVOR_SUMMARY.out.csv).join(SURVIVOR_VCF_TO_TABLE.out.annotation)
    SURVIVOR_ANNOTATION(surv_annot_input)
    surv_inexon_input = SURVIVOR_MERGE.out.vcf.join(SURVIVOR_BED_INTERSECT.out.intersected_exons)
    SURVIVOR_INEXON(surv_inexon_input)
}
