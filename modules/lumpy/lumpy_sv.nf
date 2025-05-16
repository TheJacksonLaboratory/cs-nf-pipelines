process LUMPY_SV {
    tag "$sampleID"
    
    cpus = 1
    memory 80.GB
    time '18:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/lumpy-sv:0.3.1--hdfd78af_3'
    
    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name)

    output:
    tuple val(sampleID), path("*_lumpy_sv.vcf"), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name), val('lumpy'), emit: lumpy_sv_vcf

    script:
    """
    lumpyexpress \
        -B ${tumor_bam},${normal_bam} \
        -o ${sampleID}_lumpy_sv.vcf
    """
}
