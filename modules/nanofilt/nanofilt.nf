process NANOFILT{
    tag "$sampleID" 

    cpus 16
    memory 24.GB
    time "24:00:00"
    publishDir "${params.pubdir}/", pattern: "*_porechop_NanoFilt.fastq", mode:'copy'

    container 'quay.io/biocontainers/nanofilt:2.8.0--py_0'

    input:
        tuple val(sampleID), file(porechop_fastq)

    output:
        tuple val(sampleID), file("*_porechop_NanoFilt.fastq"), emit: porechop_nanofilt_fastq

    script:
        """
        NanoFilt -q ${params.quality} -l ${params.length} --headcrop ${params.headcrop} --tailcrop ${params.tailcrop} --logfile  ${sampleID}_nonofilt_log  ${porechop_fastq}  > ${sampleID}_porechop_NanoFilt.fastq
        """
}