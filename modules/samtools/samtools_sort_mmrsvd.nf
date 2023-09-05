process SAMTOOLS_SORT {
    tag "${sampleID}"

    cpus 8
    memory 125.GB
    time '48:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/alignments' : 'alignments' }", mode:'copy', enabled: params.keep_intermediate ? true : false

    input:
        tuple val(sampleID), file(sam)

    output:
        tuple val(sampleID), file("${sampleID}.bam"), emit: bam

    script:
        """
        samtools sort --threads ${task.cpus} ${sam} -o ${sampleID}.bam
        """
}
