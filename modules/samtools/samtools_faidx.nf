process SAMTOOLS_FAIDX {
    tag "${fasta}"

    cpus 1
    memory 4.GB
    time '2:00:00'
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    input:
        file(fasta)

    output:
        tuple file("${fasta}"), file("${fasta}.fai"), emit: fasta_fai

    script:
        """
        samtools faidx ${fasta}
        """
}