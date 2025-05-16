process LOG_PARSER {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.summary_QC_metrics.txt", mode: 'copy'

    container 'python:3.8.10'

    input:
    tuple val(sampleID), file(log_cutadapt), file(log_bowtie), file(log_sorted_metrics), file(log_mtdna_content), file(log_pbc_qc), file(log_fraction_reads)

    output:
    tuple val(sampleID), file("*.summary_QC_metrics.txt")

    script:
    """
    python ${projectDir}/bin/atac/LogParser.py > ${sampleID}.summary_QC_metrics.txt
    """
}
