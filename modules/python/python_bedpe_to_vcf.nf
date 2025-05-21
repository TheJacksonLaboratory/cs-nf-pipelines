process PYTHON_BEDPE_TO_VCF {
    tag "$sampleID"

    cpus 1
    memory 20.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/fuc:0.37.0--pyh7cba7a3_0'

    publishDir "${params.pubdir}/${sampleID}", pattern:"*.vcf", mode:'copy'

    input:
        tuple val(sampleID), path(bedpe)
    output:
        tuple val(sampleID), path("*.vcf"), emit: vcf
    script:
    """
    python ${projectDir}/bin/germline_sv/bedpetovcf.py -i ${bedpe} -v ${sampleID}.mergedCall.vcf
    """
}
