process FRIP_READS_IN_PEAKS {
    tag "$sampleID"

    cpus 2
    memory 4.GB 
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
    tuple val(sampleID), file(processed_bams), file(narrow_peaks)

    output:
    tuple val(sampleID), file("reads_in_peaks.tmp.ba*")

    script:
    """
    bedtools sort \
    -i ${narrow_peaks} \
    | bedtools merge -i stdin \
    | bedtools intersect -u -nonamecheck \
    -a ${processed_bams[0]} \
    -b stdin \
    -ubam \
    > reads_in_peaks.tmp.bam
    """
}
