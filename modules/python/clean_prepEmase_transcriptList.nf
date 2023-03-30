process CLEAN_TRANSCRIPT_LISTS {

    cpus = 1
    memory = 1.GB
    time = '00:45:00'

    container 'quay.io/jaxcompsci/py3_perl_pylibs:v1'

    publishDir "${params.pubdir}/emase", pattern: '*.info', mode:'copy'

    input:
    path(pooled_transcript_list)

    output:
    path("*.info"), emit: transcript_info

    script:
    """
    python ${projectDir}/bin/emase/clean_transcript_info.py \
    --input-transcript-list ${pooled_transcript_list} \
    --haplotype-transcript-output emase.pooled.fullTranscripts.info \
    --full-transcript-output emase.fullTranscripts.info \
    --haplotype-list ${params.haplotype_list}
    """

    stub:
    """
    touch emase.pooled.fullTranscripts.info
    touch emase.fullTranscripts.info
    """

}