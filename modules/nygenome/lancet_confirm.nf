process LANCET_CONFIRM {
    tag "$sampleID"

    cpus = 8
    memory = 50.GB
    time = '20:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/lancet:v1.1.0'
    // publishDir "${params.pubdir}/${sampleID}", pattern:".vcf", mode:'copy'

    input:
    tuple val(sampleID), path(bed), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name), val(chrom)

    output:
    tuple val(sampleID), path("*.vcf"), val(meta), val(normal_name), val(tumor_name), val(chrom), emit: vcf

    script:
    """
    lancet \
    --tumor ${tumor_bam} \
    --normal ${normal_bam} \
    --bed ${bed} \
    --ref ${params.ref_fa} \
    --min-k 11 \
    --low-cov 1 \
    --min-phred-fisher 5 \
    --min-strand-bias 1 \
    --min-alt-count-tumor 3 \
    --min-vaf-tumor 0.04 \
    --padding 250 \
    --window-size 2000 \
    --num-threads ${task.cpus} \
    > ${sampleID}_lancet_merged_${chrom}.vcf
    """
}
