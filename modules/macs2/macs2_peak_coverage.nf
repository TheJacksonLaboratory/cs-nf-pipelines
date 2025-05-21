process PEAK_COVERAGE {
    tag "$sampleID"

    cpus = 1
    memory 1.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/macs2:2.2.7.1--py39hbf8eff0_4'  

    input:
    tuple val(sampleID), file(narrow_peaks)

    output:
    tuple val(sampleID), file("*_peaks.narrowPeak.saf")

    shell:
    '''
    awk 'OFS="\\t" {print $1"."$2"."$3, $1, $2, $3, "."}' !{narrow_peaks} \
    > !{sampleID}_peaks.narrowPeak.saf
    '''
}
