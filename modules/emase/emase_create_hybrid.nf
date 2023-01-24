process EMASE_CREATE_HYBRID {
    
    // give a group of fastas, and a haplotype list 
    // 1. generate a hybrid genome
    // 2. generate transcript list
    // 3. generate bowtie index. 

    cpus 1
    memory {15.GB * task.attempt}
    time {24.hour * task.attempt}
    errorStrategy 'retry' 
    maxRetries 1

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:3ac8573'

    publishDir "${params.pubdir}/emase", pattern: "*.fa", mode:'copy'
    publishDir "${params.pubdir}/emase", pattern: "*.info", mode:'copy'
    publishDir "${params.pubdir}/emase", pattern: "*.tsv", mode:'copy'
    publishDir "${params.pubdir}/emase/bowtie", pattern: "*.ebwt", mode:'copy'

    output:
    path file("*.fa"), emit: transcript_fasta
    path file("*.info"), emit: transcript_info
    path file("*.ebwt"), emit: bowtie_index

    script:
    """
    create-hybrid -F ${params.genome_file_list} -s ${params.haplotype_list} -o ./ --create-bowtie-index
    """

    stub:
    """
    touch emase.pooled.transcripts.fa
    touch emase.pooled.transcripts.info
    touch emase.gene2transcripts.tsv
    touch bowtie.transcripts.4.ebwt
    touch bowtie.transcripts.3.ebwt
    touch bowtie.transcripts.2.ebwt
    touch bowtie.transcripts.1.ebwt
    touch bowtie.transcripts.rev.2.ebwt
    touch bowtie.transcripts.rev.1.ebwt
    """
}


/*
create-hybrid:
Usage:
    create-hybrid -F <fasta_files> [ -s <hap_list> -o <out_file> --create-bowtie-index ]

Input:
    -F <fasta_files> : List of fasta files (comma delimited)
    -s <hap_list>    : Names of haplotypes to be used instead (comma delimited, in the order of genomes)
    -o <out_file>    : Output file name (default: './emase.pooled.targets.fa')

Parameters:
    --help, -h            : shows this help message
    --create-bowtie-index : builds bowtie1 index
*/