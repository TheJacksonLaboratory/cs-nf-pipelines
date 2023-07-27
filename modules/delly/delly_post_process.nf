process DELLY_POST_PROCESS {
    tag "$sampleID"
    
    cpus = 1
    memory = 20.GB
    time = "2:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.10.2--h4f4756c_3'
    
    input:
        tuple val(sampleID), file(delly_bcf)
    
    output:
        tuple val(sampleID), file("${sampleID}_delly.vcf"), emit: delly_vcf

    script:
        """
        bcftools view ${delly_bcf} > "${sampleID}_delly_unsorted.vcf"
        bash ${projectDir}/bin/vcfSort.sh "${sampleID}_delly_unsorted.vcf" ${sampleID}_delly.vcf
        """
}