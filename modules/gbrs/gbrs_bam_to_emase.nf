process BAM_TO_EMASE {
    tag "$sampleID"

    cpus 8
    memory {60.GB * task.attempt}
    time {30.hour * task.attempt}
    errorStrategy 'retry' 
    maxRetries 1

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:89bbb10'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'bowtie' }", pattern: "*.h5", mode: 'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), file("*.emase.h5"), emit: emase_h5

    script:
    """
    gbrs bam2emase -i ${bam} \
                -m ${params.transcripts_info} \
                -s ${params.gbrs_strain_list} \
                -o ${bam.baseName}_emase.h5
    """

    stub:
    """
    touch ${bam.baseName}_emase.h5
    """
}


