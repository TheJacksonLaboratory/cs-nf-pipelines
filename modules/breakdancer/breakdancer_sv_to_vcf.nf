process BREAKDANCER_SV_TO_VCF {
    tag "$sampleID"
    
    cpus = 1
    memory = 20.GB
    time = "2:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'
    
    input:
        tuple val(sampleID), file(bam), file(bai), file(breakdancer_sv)
    
    output:
        tuple val(sampleID), file("${sampleID}_breakdancer_sorted.vcf"), emit: breakdancer_vcf

    script:
        """
        python ${projectDir}/bin/breakdancer2vcfHeader.py -i ${breakdancer_sv} -o ${sampleID}_breakdancer.vcf
        bash ${projectDir}/bin/vcfSort.sh ${sampleID}_breakdancer.vcf ${sampleID}_breakdancer_sorted.vcf
        """
}