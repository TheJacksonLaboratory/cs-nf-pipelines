process CLEAN_TRANSCRIPT_LISTS {

    cpus = 1
    memory = 1.GB
    time = '00:45:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

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
