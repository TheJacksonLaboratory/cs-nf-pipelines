process GZIP {
    tag "$sampleID"

    cpus 1  
    memory 15.GB
    time 18.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/jaxcompsci/py3_perl_pylibs:v2"

    publishDir "${params.pubdir}/${sampleID + '/processed_reads'}", pattern: "*.gz", mode:'copy'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("*.gz"), emit: gunzip_fastq
    
    script:
    """
    gzip -c ${reads[0]} > ${reads[0]}.gz
    gzip -c ${reads[1]} > ${reads[1]}.gz
    """
}
