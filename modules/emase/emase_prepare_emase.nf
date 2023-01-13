process EMASE_PREPARE_EMASE {

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

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:89bbb10'

    publishDir "${params.pubdir}/emase", pattern: "[*.fa, *.info, *.tsv]", mode:'copy'
    publishDir "${params.pubdir}/emase/bowtie", pattern: "[*.ebwt]", mode:'copy'

    output:
    path file("*.fa"), emit: transcript_fasta
    path file("*.info"), emit: transcript_info
    path file("*.tsv"), emit: gene_to_transcripts
    path file("*.ebwt"), emit: bowtie_index

    script:
    """
    prepare-emase -G ${params.genome_file_list} -g ${params.gtf_file_list} -s ${params.haplotype_list} -o ./ -m
    """
}