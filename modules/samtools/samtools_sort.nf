process SAMTOOLS_SORT {
    tag "$sampleID"

    cpus 4
    memory 20.GB
    time '20:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*.bam", mode:'copy', enabled: params.workflow == 'rrbs' ? true : false
    
    input:
    tuple val(sampleID), file(sam_file)
    val(options)
    val(suffix)

    output:
    tuple val(sampleID), file("*.sorted.*"), emit: sorted_file

    script:
    """
    samtools sort \
    ${options} \
    -@ ${task.cpus} \
    -o ${sam_file.baseName}.sorted.${suffix} \
    ${sam_file}
    """
}
