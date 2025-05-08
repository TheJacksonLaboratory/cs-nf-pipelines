process NANOFILT{
    tag "$sampleID" 

    cpus 16
    memory 24.GB
    time "24:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    publishDir "${params.pubdir}/${sampleID + '/fastq'}", pattern: "*_porechop_NanoFilt.fastq", mode:'copy', enabled: params.keep_intermediate ? true : false

    container 'quay.io/biocontainers/nanofilt:2.8.0--py_0'

    input:
        tuple val(sampleID), file(porechop_fastq)

    output:
        tuple val(sampleID), file("*_porechop_NanoFilt.fastq"), emit: porechop_nanofilt_fastq

    script:
        """
        NanoFilt -q ${params.quality} -l ${params.length} --headcrop ${params.headcrop} --tailcrop ${params.tailcrop} --logfile  ${sampleID}_nonofilt_log  ${porechop_fastq}  > ${sampleID}_porechop_NanoFilt.fastq
        """
}
