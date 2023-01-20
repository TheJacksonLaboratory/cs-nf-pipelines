process ALNTOOLS_BAM2EMASE {
    tag "$sampleID"

    cpus 1
    memory {15.GB * task.attempt}
    time {10.hour * task.attempt}
    errorStrategy 'retry' 
    maxRetries 1

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:3ac8573'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'alntools' }", pattern: "*.h5", mode: 'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), file("*.emase.h5"), emit: emase_h5

    script:
    """
    alntools bam2emase \
                -t ${params.transcripts_info} \
                ${bam} \
                ${bam.baseName}.emase.h5
    """

    stub:
    """
    touch ${bam.baseName}.emase.h5
    """
}


