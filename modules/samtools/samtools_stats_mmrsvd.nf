process SAMTOOLS_STATS {
    tag "${sampleID}"

    cpus 1
    memory 2.GB
    time '1:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    publishDir "${params.pubdir}/${sampleID + '/alignments'}", mode:'copy'

    input:
        tuple val(sampleID), file(bam), file(bai)

    output:
        tuple val(sampleID), file("${sampleID}.insert_size.txt"), emit: insert_size

    script:
        """
        samtools stats ${bam} |grep "^IS" |awk '{a = a + \$2*\$3; b = b + \$3}END{print int(a/b)}' > ${sampleID}.insert_size.txt
        """
}
