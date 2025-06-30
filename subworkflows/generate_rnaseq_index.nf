#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/generate_rnaseq_index.nf"
include {param_log} from "${projectDir}/bin/log/generate_rnaseq_index.nf"
include {MAKE_CUSTOM_TRANSCRIPTOME} from "${projectDir}/modules/utility_modules/make_custom_transcriptome"
include {RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_BOWTIE2;
         RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_STAR} from "${projectDir}/modules/rsem/rsem_preparereference"
include {PICARD_CREATESEQUENCEDICTIONARY} from "${projectDir}/modules/picard/picard_createsequencedictionary"
include {UCSC_GTFTOGENEPRED} from "${projectDir}/modules/ucsc/ucsc_gtftogenepred"
include {GENERATE_RRNA_INTERVALS} from "${projectDir}/modules/utility_modules/generate_rrna_intervals"
include {KALLISTO_INDEX} from "${projectDir}/modules/kallisto/kallisto_index"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

def checkFileExists(filePath, name) {
    if (filePath && !file(filePath).exists()) {
        log.error "File not found: ${filePath} (${name})"
        exit 1
    }
}

checkFileExists(params.ref_fa, "ref_fa")
checkFileExists(params.ref_gtf, "ref_gtf")
if (params.custom_gene_fasta) {
    checkFileExists(params.custom_gene_fasta, "custom_gene_fasta")
}

workflow GENERATE_RNASEQ_INDEX {

    if (params.custom_gene_fasta) {
    
        MAKE_CUSTOM_TRANSCRIPTOME(Channel.of([params.ref_fa, params.ref_gtf, params.custom_gene_fasta]))

        bowtie2_input = MAKE_CUSTOM_TRANSCRIPTOME.out.concat_fasta.combine(MAKE_CUSTOM_TRANSCRIPTOME.out.concat_gtf).map{it -> [it[0], it[1], 'bowtie2', '']}
        star_build_set = MAKE_CUSTOM_TRANSCRIPTOME.out.concat_fasta.combine(MAKE_CUSTOM_TRANSCRIPTOME.out.concat_gtf).combine(Channel.of(75, 100, 125, 150)).map{it -> [it[0], it[1], 'STAR', it[2]]}
        
        fasta = MAKE_CUSTOM_TRANSCRIPTOME.out.concat_fasta
        gtf = MAKE_CUSTOM_TRANSCRIPTOME.out.concat_gtf
    
    } else {
        bowtie2_input = Channel.of([params.ref_fa, params.ref_gtf, 'bowtie2', ''])
        star_build_set = Channel.of([params.ref_fa, params.ref_gtf, 'STAR']).combine(Channel.of(75, 100, 125, 150))
        
        fasta = params.ref_fa
        gtf = params.ref_gtf
    }

    RSEM_PREPAREREFERENCE_BOWTIE2(bowtie2_input)
    
    RSEM_PREPAREREFERENCE_STAR(star_build_set)

    PICARD_CREATESEQUENCEDICTIONARY(fasta)

    GENERATE_RRNA_INTERVALS(gtf, PICARD_CREATESEQUENCEDICTIONARY.out.dict)
    
    UCSC_GTFTOGENEPRED(gtf)

    KALLISTO_INDEX(RSEM_PREPAREREFERENCE_BOWTIE2.out.transcripts)

}

// Workflow adapted from: https://github.com/KU-GDSC/workflows
