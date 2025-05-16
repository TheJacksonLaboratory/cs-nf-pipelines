process BAM_COVERAGE_BIGWIG {
    tag "$sampleID"

    cpus 8
    memory 10.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/deeptools'}", pattern: "*.bigwig", mode: 'copy'
    container 'quay.io/biocontainers/deeptools:3.3.2--py_1'

    input:
    tuple val(sampleID), file(processed_bams)

    output:
    tuple val(sampleID), file("*.bigwig")

    script:
    """
    bamCoverage \
    --numberOfProcessors $task.cpus \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize ${params.effective_genome_size} \
    --bam ${processed_bams[0]} \
    --outFileFormat bigwig \
    --outFileName ${sampleID}.bigwig
    """
}
