process REHEADER_VCF {
    tag "$sampleID"

    cpus = 1
    memory = 1.GB
    time '00:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_3"

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/unmerged_calls' : 'unmerged_calls' }", pattern: "${sampleID}_${caller}Sort.vcf", mode: 'copy'

    input:
        tuple val(sampleID), file(vcf)
        val(caller)

    output:
        tuple val(sampleID), file("${sampleID}_${caller}Sort.vcf"), emit: vcf_rehead

    script:
        """
        printf "${sampleID}_${caller}\n" > rehead.txt
        bcftools reheader --samples rehead.txt \
            -o ${sampleID}_${caller}Sort.vcf \
            ${vcf}
        """	
}