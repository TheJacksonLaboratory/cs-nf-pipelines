process LUMPY_COMPRESS_INDEX {
    tag "$sampleID"

    cpus = 1
    memory = 2.GB
    time = '00:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'
    
    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern:"*.gz", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), val(meta), val(normal_name), val(tumor_name), val(caller)

    output:
    tuple val(sampleID), path("*sorted.vcf.gz"), path("*.tbi"), val(meta), val(normal_name), val(tumor_name), val(caller), emit: vcf_tbi

    script:
    """
    cat ${vcf} | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > ${vcf.baseName}_sorted.vcf
    bgzip ${vcf.baseName}_sorted.vcf
    tabix -p vcf ${vcf.baseName}_sorted.vcf.gz
    """
}
