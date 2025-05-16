process BCFTOOLS_MERGE_DELLY_CNV {
    tag "$sampleID"
    
    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    tuple val(sampleID), path(normal_bcf), path(tumor_bcf), val(meta), val(normal_name), val(tumor_name)

    output:
    tuple val(sampleID), path("${sampleID}.bcf"), path("${sampleID}.bcf.csi"), val(meta), val(normal_name), val(tumor_name), emit: merged_bcf

    script:

    """
    bcftools index ${tumor_bcf}
    bcftools index ${normal_bcf}
    bcftools merge -m id -O b -o ${sampleID}.bcf ${tumor_bcf} ${normal_bcf}
    bcftools index ${sampleID}.bcf
    """
}
