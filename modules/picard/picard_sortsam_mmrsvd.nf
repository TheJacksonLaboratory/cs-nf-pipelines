process PICARD_SORTSAM {
    tag "${sampleID}"

    cpus 1
    memory {bam.size() < 40.GB ? 20.GB : 40.GB}
    time {bam.size() < 40.GB ? '3:00:00' : '12:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/picard:2.25.1--hdfd78af_1'

    input:
        tuple val(sampleID), file(bam)
        val(suffix_string)

    output:
        tuple val(sampleID), file("${sampleID}_${suffix_string}.bam"), file("${sampleID}_${suffix_string}.bai"), emit: sorted_bam


    script:
        """
        picard SortSam I=${bam} O=${sampleID}_${suffix_string}.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
        """
}
