process NANOQC{
    tag "$sampleID" 

    cpus 16
    memory 24.GB
    time "24:00:00"

    publishDir "${params.pubdir}/stats", pattern: "*_porechop_1st_nanoQC", mode:'copy'

    container 'quay.io/biocontainers/nanoqc:0.9.4--py_0'

    input:
        tuple val(sampleID), file(porechop_fastq)

    output:
        tuple val(sampleID), file("*_porechop_1st_nanoQC"), emit: porechop_nanoqc

    script:
        """
        nanoQC -o ${sampleID}_porechop_1st_nanoQC ${porechop_fastq}
        """
}