process BUILD_BWA_INDEX {
    tag "${sampleID}"

    cpus 1
    memory 30.GB
    time '02:00:00'
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_6'
    publishDir ${params.outdir}, mode: 'copy'

    input:
        file(fasta)
    output:
        file("${fasta}.*"), emit: bwa_index

    script:
        """
        bwa index ${fasta}
        """
}