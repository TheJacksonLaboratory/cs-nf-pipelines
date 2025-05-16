process SAMTOOLS_FILTER {
    tag "$sampleID"

    cpus 2
    memory 4.GB
    time '10:00:00'

    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    input:
    tuple val(sampleID), file(in_file)
    val(option)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam

    script:
    // Exclude reads based on input bit flag.

    prefix = "${sampleID}.Lb"
    if (params.workflow == "chipseq"){
        output = "${prefix}.bam"
    }
    else{
        output = "${sampleID}.bam"
    }
    """
    samtools view -h -b ${option} ${in_file} > ${output}
    """
}
