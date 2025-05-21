process NANOSTAT {
    tag "$sampleID"

    cpus 16
    memory 24.GB
    time "24:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "nanostat*", mode:'copy'

    container 'quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0'

    input:
        tuple val(sampleID), file(read1)

    output:
        path("nanostat_*"), emit: nanostat_report

    script:
        """
        NanoStat --fastq ${read1} -n nanostat_${read1.baseName}_${sampleID}
        """
}
