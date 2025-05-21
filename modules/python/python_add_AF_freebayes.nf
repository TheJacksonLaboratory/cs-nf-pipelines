process ADD_ALT_AF {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/bedtools-python3:2.26.0'

    publishDir "${params.pubdir}/${sampleID}", pattern:"*.vcf", mode:'copy'

    input:
    tuple val(sampleID), path(vcf)
    val(output_suffix)

    output:
    tuple val(sampleID), path("*.vcf"), emit: vcf

    script:
    """
    python ${projectDir}/bin/wes/AF_freebayes.py ${vcf} ${sampleID}_${output_suffix}.vcf 9
    """
}
