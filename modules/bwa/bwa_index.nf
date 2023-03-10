process BWA_INDEX {
    tag "${fasta}"

    cpus 1
    memory 30.GB
    time '02:00:00'
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_6'

    input:
        file(fasta)
    output:
        path("${fasta}.*"), emit: bwa_index

    script:
        """
        bwa index ${fasta}
        """
}