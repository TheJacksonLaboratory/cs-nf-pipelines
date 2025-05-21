process BAMTOOLS_STATS {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/bamtools:2.5.1--h9a82719_9'

    publishDir "${params.pubdir}/${sampleID}", pattern:"*.txt", mode:'copy'

    input:
    tuple val(sampleID), file(reordered_sorted_bam)

    output:
    tuple val(sampleID), file("*metrics.txt"), emit: picard_metrics

    script:
    if (params.read_type == "PE")

        """
        bamtools stats -insert -in ${reordered_sorted_bam} > ${sampleID}_aln_metrics.txt
        """

    else if (params.read_type == "SE")

        """
        bamtools stats -in ${reordered_sorted_bam} > ${sampleID}_aln_metrics.txt
        """
}
