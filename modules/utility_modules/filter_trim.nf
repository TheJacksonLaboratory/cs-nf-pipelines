process FILTER_TRIM {
    tag "$sampleID"

    cpus 1
    memory 30.GB
    time '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

    publishDir "${params.pubdir}/stats", pattern: "*_stat", mode:'copy'

    input:
    tuple val(sampleID), file(fq_reads)

    output:
    tuple val(sampleID), file("*_stat"), emit: quality_stats
    tuple val(sampleID), file("*filtered_trimmed"), emit: trimmed_fastq

    script:

    if (!params.fastq2){
        inputfq="-S ${fq_reads[0]}"
    }
    else {
        inputfq="${fq_reads[0]} ${fq_reads[1]}"
    }

    """
    python ${projectDir}/bin/germline_sv/filter_trim.py $inputfq
    """
}
