#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/joint_gvcf_calling.nf"
include {param_log} from "${projectDir}/bin/log/joint_gvcf_calling.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv_gvcf.nf"
include {GATK_INDEXFEATUREFILE} from "${projectDir}/modules/gatk/gatk_indexfeaturefile"
include {GATK_COMBINEGVCFS_INTERVALS} from "${projectDir}/modules/gatk/gatk_combinegvcfs_intervals"
include {GATK_GENOTYPEGVCF} from "${projectDir}/modules/gatk/gatk_genotypegvcf"
include {GATK_GATHERVCFS} from "${projectDir}/modules/gatk/gatk_gathervcfs"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// prepare reads channel
if (params.csv_input) {
    def csvFile = file(params.csv_input)
    if (!csvFile.exists()) {
        exit 1, "ERROR: CSV input file does not exist. Parameter csv_input was set to: ${params.csv_input}"
    }
    ch_input_sample = extract_csv(csvFile)
    ch_input_sample.map{it -> [it[0], it[2]]}.set{gvcf_ch}
    ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
} else {
    exit 1, "ERROR: CSV input is required. Parameter csv_input was set to: ${params.csv_input}"
}

// main workflow
workflow JOINT_GVCF_CALLING {

    GATK_INDEXFEATUREFILE(gvcf_ch)

    // Read a list of contigs from parameters to provide to GATK as intervals
    chroms = Channel
        .fromPath("${params.chrom_contigs}")
        .splitText()
        .map{it -> it.trim()}
    num_chroms = file(params.chrom_contigs).countLines().toInteger()

    combine_input_ch = gvcf_ch
    .join(GATK_INDEXFEATUREFILE.out.idx) 
    .join(meta_ch)
    .map{ it -> [it[3]['output'], it[1], it[2]] }
    .groupTuple() // make this sized by number of samples from meta_ch
    .combine(chroms)

    GATK_COMBINEGVCFS_INTERVALS(combine_input_ch)

    GATK_GENOTYPEGVCF(GATK_COMBINEGVCFS_INTERVALS.out.gvcf_idx)

    GATK_GATHERVCFS(GATK_GENOTYPEGVCF.out.vcf.groupTuple(size: num_chroms), 'Genotyped_allChroms')
}
