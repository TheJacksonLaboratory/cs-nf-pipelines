process RENAME_VCF {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/bedtools-python3:2.26.0'

    input:
    tuple val(sampleID), file(vcf), val(meta), val(normal_name), val(tumor_name), val(tool)

    output:
    tuple val(sampleID), file("*_sampleNamed.vcf"), val(meta), val(normal_name), val(tumor_name), val(tool), emit: rename_vcf

    script:
    
    normal = meta.normal_id
    tumor = meta.tumor_id

    tool_name = tool == 'lancet_support' ? 'lancet' : tool

    """
    python \
    ${projectDir}/bin/pta/rename_vcf.py \
    ${vcf} \
    ${vcf.baseName}_sampleNamed.vcf \
    ${normal} \
    ${tumor} \
    ${tool} 
    """
}
