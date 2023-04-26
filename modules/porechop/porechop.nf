process PORECHOP {
    tag "$sampleID"

    cpus 16
    memory 24.GB
    time "24:00:00"

    publishDir "${params.pubdir}/fastq", pattern: "${sampleID}_porechop.fastq", mode:'copy', enabled: params.keep_intermediate ? true : false

    container 'quay.io/biocontainers/porechop:0.2.4--py39hc16433a_3'

    input:
        tuple val(sampleID),file(read1)

    output:
        tuple val(sampleID), file("${sampleID}_porechop.fastq"), emit: porechop_fastq

    script:
        """
        porechop -i ${read1} -o ${sampleID}_porechop.fastq  --format fastq -t ${task.cpus}
        """
}