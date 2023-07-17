process BOWTIE_BUILD {
    cpus 8
    memory {30.GB}
    time {15.hour}

    container 'quay.io/jaxcompsci/bowtie-samtools:v1'

    publishDir "${params.pubdir}/bowtie", pattern: '*.ebwt', mode:'copy'

    input:
    path(fasta)
    val(output_prefix)

    output:
    path("*.ebwt"), emit: bowtie_ref

    script:
    """
    bowtie-build --threads ${task.cpus} ${fasta} ${output_prefix}
    """

    stub:
    """
    touch "bowtie.transcripts.1.ebwt"
    touch "bowtie.transcripts.2.ebwt"
    touch "bowtie.transcripts.3.ebwt"
    touch "bowtie.transcripts.4.ebwt"
    """
}
