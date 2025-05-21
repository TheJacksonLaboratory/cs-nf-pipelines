process BISMARK_DEDUPLICATION {
    tag "$sampleID"

    cpus 8
    memory 60.GB
    time 30.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bismark:0.23.1--hdfd78af_0'

    publishDir "${params.pubdir}/${sampleID + '/alignment'}", pattern: "*.bam", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*txt", mode:'copy'

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*.bam"), emit: dedup_bam
    tuple val(sampleID), file("*report.txt"), emit: dedup_report

    script:

    fq_type = params.read_type == 'PE' ? '-p' : '-s'
    
    """
    deduplicate_bismark $fq_type --bam $bam
    """
}
