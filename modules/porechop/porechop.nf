process PORECHOP {
    tag "$sampleID"

    cpus 20
    memory 200.GB
    time "72:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/fastq'}", pattern: "${sampleID}_porechop.fastq", mode:'copy', enabled: params.keep_intermediate ? true : false

    container 'quay.io/biocontainers/porechop:0.2.4--py39hc16433a_3'

    input:
        tuple val(sampleID),file(read1)

    output:
        tuple val(sampleID), file("${sampleID}_porechop.fastq"), emit: porechop_fastq

    script:
        """
        porechop -i ${read1} -o ${sampleID}_porechop.fastq  --format fastq -t ${task.cpus}
        """
}
