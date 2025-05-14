process LYMPHOMA_CLASSIFIER {
    tag "$sampleID"

    cpus 1
    memory 1.GB
    time '00:20:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'

    publishDir "${params.pubdir}/${sampleID}", pattern:"*EBV_classifier.txt", mode:'copy'

    input:
    tuple val(sampleID), file(rsem_counts)

    output:
    tuple val(sampleID), file("*EBV_classifier.txt"), emit: ebv_classification

    script:
    """
    python \
    ${projectDir}/bin/rnaseq/lymphoma_classifier.py \
    --expected_expression ${params.classifier_table} \
    --counts ${rsem_counts} \
    --output ${sampleID}.EBV_classifier.txt
    """
}