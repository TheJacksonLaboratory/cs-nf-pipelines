process BCFTOOLS_COMPRESS_INDEX {
    tag "$sampleID"
    
    cpus = 1
    memory = 2.GB
    time = '00:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'
    
    publishDir "${params.pubdir}/${sampleID}", pattern: "*haplotypecaller.gatk.filtered.vcf.gz", mode:'copy'

    input:
    tuple val(sampleID), file(vcf)

    output:
    tuple val(sampleID), file("*.gz"), file("*.gz.tbi"), emit: vcf_idx

    script:
    
    name_adjust = params.gen_org == 'mouse' && params.workflow =='pta' ? "mv ${vcf}.gz ${sampleID}_haplotypecaller.gatk.filtered.vcf.gz; mv ${vcf}.gz.tbi ${sampleID}_haplotypecaller.gatk.filtered.vcf.gz.tbi" : ""

    """
    bgzip -f ${vcf}
    tabix -p vcf ${vcf}.gz
    ${name_adjust}
    """
}
