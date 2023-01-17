process EMASE_PREPARE_TRANSCRIPT_LIST {

    // give a fasta or group of fastas, gtf or group of gtfs, and a haplotype list 
    // 1. generate a hybrid genome
    // 2. generate transcript list
    // 3. generate gene to transcript map
    // 4. generate bowtie index. 

    cpus 1
    memory {60.GB * task.attempt}
    time {30.hour * task.attempt}
    errorStrategy 'retry' 
    maxRetries 1

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:daafe97'

    publishDir "${params.pubdir}/emase", pattern: '*.fa', mode:'copy'
    publishDir "${params.pubdir}/emase", pattern: '*.info', mode:'copy'
    publishDir "${params.pubdir}/emase", pattern: '*.tsv', mode:'copy'

    output:
    path("*.fa"), emit: transcript_fasta
    path("*.info"), emit: transcript_info
    path("*.tsv"), emit: gene_to_transcripts

    script:
    """
    prepare-emase -G ${params.base_genome} -g ${params.base_gtf} -o ./ --no-bowtie-index
    """

    stub:
    """
    touch emase.transcripts.fa
    touch emase.transcripts.info
    touch emase.gene2transcripts.tsv
    """

}