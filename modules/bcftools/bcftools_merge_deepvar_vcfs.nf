process BCFTOOLS_MERGEDEEPVAR {
    tag "$sampleID"

    cpus 2
    memory 50.GB
    time {5.hour * task.attempt}
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'retry'}
    maxRetries 2

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.*vcf*", mode:'copy'

    input:
    tuple val(sampleID), path(files), path(index)
    val(suffix)

    output:
    tuple val(sampleID), path("*_sorted_deepvariant.*.gz"), path("*_sorted_deepvariant.*.gz.tbi"), emit: vcf_idx

    script:
    """
    ## GVCF
    # concatenate
    bcftools concat ${files} --allow-overlaps --remove-duplicates -Oz -o ${sampleID}_deepvariant.${suffix}.gz 

    # sort
    bcftools sort ${sampleID}_deepvariant.${suffix}.gz -Oz -o ${sampleID}_sorted_deepvariant.${suffix}.gz

    # index
    bcftools index --tbi ${sampleID}_sorted_deepvariant.${suffix}.gz
    """
}
