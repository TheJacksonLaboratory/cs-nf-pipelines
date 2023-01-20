process EMASE_GET_COMMON_ALIGNMENTS {

    cpus 1
    memory {15.GB * task.attempt}
    time {5.hour * task.attempt}
    errorStrategy 'retry' 
    maxRetries 1

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:3ac8573'

    publishDir "${params.pubdir}/emase", pattern: '*.fa', mode:'copy'
    publishDir "${params.pubdir}/emase", pattern: '*.info', mode:'copy'
    publishDir "${params.pubdir}/emase", pattern: '*.tsv', mode:'copy'
    publishDir "${params.pubdir}/emase/bowtie", pattern: "*.ebwt", mode:'copy'

    output:
    path("*.fa"), emit: pooled_transcript_fasta
    path("*.info"), emit: pooled_transcript_info
    path("*.tsv"), emit: pooled_gene_to_transcripts
    path("*.ebwt"), emit: pooled_bowtie_index

    script:
    """
    get-common-alignments -i ${EMASE_FILE_R1},${EMASE_FILE_R2} -o ${EMASE_FILE}
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

