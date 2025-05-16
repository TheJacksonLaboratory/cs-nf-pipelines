process LANCET {
    tag "$sampleID"

    cpus = 4
    memory = 30.GB
    time = '18:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/lancet:v1.1.0'
    publishDir "${params.pubdir}/${sampleID}", pattern:"*.vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name), path(bed), val(index)

    output:
    tuple val(sampleID), path("*_lancet.vcf"), val(meta), val(normal_name), val(tumor_name), val('lancet'), emit: vcf

    script:
    """
    lancet \
    --tumor ${tumor_bam} \
    --normal ${normal_bam} \
    --ref ${params.ref_fa} \
    --bed ${bed} \
    --min-k 11 \
    --low-cov 1 \
    --min-phred-fisher 5 \
    --min-strand-bias 1 \
    --min-alt-count-tumor 3 \
    --min-vaf-tumor 0.04 \
    --num-threads ${task.cpus} > ${sampleID}_${index}_lancet.vcf
    """
}
