process RENAME_METADATA {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/bedtools-python3:2.26.0'

    input:
    tuple val(sampleID), file(vcf), path(idx), val(meta), val(normal_name), val(tumor_name), val(tool)

    output:
    tuple val(sampleID), file("*headerAdjust.vcf"), val(meta), val(normal_name), val(tumor_name), val(tool), emit: rename_metadata_vcf

    script:
    output_name = vcf.getBaseName().replace('.vcf', '')
    """
    gunzip -c ${vcf} > temp.vcf
    python \
    ${projectDir}/bin/pta/rename_metadata.py \
    temp.vcf \
    ${output_name}_headerAdjust.vcf \
    ${tool} 
    """
}
