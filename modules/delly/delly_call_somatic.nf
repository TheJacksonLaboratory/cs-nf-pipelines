process DELLY_CALL_SOMATIC {
    tag "$sampleID"
    
    cpus = 1
    memory 80.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/delly:1.1.6--h6b1aa3f_2'
    
    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name)

    output:
    tuple val(sampleID), path("*.bcf"), path("*.csi"), val(meta), val(normal_name), val(tumor_name), emit: bcf_csi

    script:
    """
    delly call -x ${params.delly_exclusion} -o ${sampleID}_delly_somaticSV.bcf -g ${params.ref_fa} ${tumor_bam} ${normal_bam} 
    """
}
