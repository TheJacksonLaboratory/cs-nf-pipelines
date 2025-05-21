process COVCALC_GATK {
    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.bed", mode:'copy'

    input:
    tuple val(sampleID), file(txt)
    val(filename)

    output:
    tuple val(sampleID), file("*.bed"), emit: bed

    script:
    """
    python ${projectDir}/bin/rnaseq/coveragecalculator.py ${txt} ${sampleID}_${filename}_avg_median_coverage.bed
    """
}
