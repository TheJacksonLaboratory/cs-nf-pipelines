process SVTYPER {
    tag "$sampleID"

    cpus = 1
    memory 24.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'jmonlong/svtyper:release-0.7.1'

    input:
    tuple val(sampleID), path(vcf), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name), val(caller)

    output:
    tuple val(sampleID), path("*_lumpy_sv_typed.vcf"), val(meta), val(normal_name), val(tumor_name), val(caller), emit: typed_vcf

    script:
    """
    svtyper \
    -i ${vcf} \
    -B ${tumor_bam},${normal_bam} \
    -l ${sampleID}_SVTYPER_bamInfo.json \
    > ${sampleID}_lumpy_sv_typed.vcf

    """
}
