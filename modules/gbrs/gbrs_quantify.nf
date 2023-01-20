process GBRS_QUANTIFY {
    tag "$sampleID"

    cpus 1
    memory {15.GB * task.attempt}
    time {5.hour * task.attempt}
    // errorStrategy 'retry' 
    // maxRetries 1

    container 'quay.io/jaxcompsci/emase_gbrs_alntools:3ac8573'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gbrs' }", pattern: "*.h5", mode: 'copy'

    input:
    tuple val(sampleID), path(bam)

    output:
    // tuple val(sampleID), file("*.compressed.emase.h5"), emit: compressed_emase_h5

    script:
    """
    gbrs quantify \
        -i ${bam} \
        -g ${params.gene2transcript_list} \
        -L ${params.full_transcript_list} \
        -M ${params.gbrs_model} \
        --report-alignment-counts \
        -o ${sampleID}
    """

    stub:
    """
    touch ${sampleID}
    """
}