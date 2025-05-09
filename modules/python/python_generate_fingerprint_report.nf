process GENERATE_FINGERPRINT_REPORT {
    tag "$sampleID"

    cpus 1
    memory 100.MB
    time '00:20:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'

    publishDir "${params.pubdir}/${sampleID}", pattern:"*.tsv", mode:'copy'

    input:
    tuple val(sampleID), file(vcf)

    output:
    tuple val(sampleID), file("*.tsv"), emit: report

    script:
    """
    python \
    ${projectDir}/bin/amplicon/generate_fingerprint_report.py \
    --input_file ${vcf} \
    --output_prefix ${sampleID}_fingerprint \
    --rsid_file ${params.amplicon_rsid_targets}
    """
}
