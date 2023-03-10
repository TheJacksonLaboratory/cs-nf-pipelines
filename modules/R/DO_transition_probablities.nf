process DO_TRANSITION_PROBABILITIES {

    cpus 1
    memory 40.GB
    time 48.hour

    container 'quay.io/jaxcompsci/r-qtl2-deseq-biomart-tidy:v1'

    // publishDir "${params.pubdir}/emase", pattern: '*.fa', mode:'copy'
    // publishDir "${params.pubdir}/emase", pattern: '*.info', mode:'copy'
    // publishDir "${params.pubdir}/emase", pattern: '*.tsv', mode:'copy'
    // publishDir "${params.pubdir}/emase/bowtie", pattern: "*.ebwt", mode:'copy'

    // output:
    // path("*.fa"), emit: pooled_transcript_fasta
    // path("*.info"), emit: pooled_transcript_info
    // path("*.tsv"), emit: pooled_gene_to_transcripts
    // path("*.ebwt"), emit: pooled_bowtie_index

    script:
    """
    Rscript ${projectDir}/bin/gbrs/gene_bp_to_cM_to_transprob.R
    """

    stub:
    """

    """

}

/*

*/