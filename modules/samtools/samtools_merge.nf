process SAMTOOLS_MERGE {
    tag "$sampleID"

    cpus 8
    memory 60.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    publishDir "${params.pubdir}/${sampleID}", pattern:"*.bam", mode:'copy', enabled: params.keep_intermediate
    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern:"*.bam", mode:'copy', enabled: params.merge_inds

    input:
        tuple val(sampleID), file(bam)
        val(filename)

    output:
        tuple val(sampleID), file("*.bam"), emit: bam

    script:
        """
        samtools merge -@ ${task.cpus} ${sampleID}_${filename}.bam ${bam} 
        """
}
