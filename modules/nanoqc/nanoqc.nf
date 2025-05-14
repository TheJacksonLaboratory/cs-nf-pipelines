process NANOQC{
    tag "$sampleID" 

    cpus 16
    memory 24.GB
    time "24:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*_porechop_1st_nanoQC", mode:'copy'

    container 'quay.io/biocontainers/nanoqc:0.9.4--py_0'

    input:
        tuple val(sampleID), file(porechop_fastq)

    output:
        tuple val(sampleID), file("*_porechop_1st_nanoQC"), emit: porechop_nanoqc

    script:
        """
        nanoQC -o ${sampleID}_porechop_1st_nanoQC ${porechop_fastq}
        """
}
