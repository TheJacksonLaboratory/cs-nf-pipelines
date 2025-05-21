
process FASTQC {
    tag "$sampleID"

    cpus 8
    memory 4.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? 'fastqc/' : '' 
        "${params.pubdir}/${type + sampleID + '/stats'}"
    }, pattern: "*_fastqc.{zip,html}", mode: 'copy'


    input:
    tuple val(sampleID), path(fq_reads)

    output:
    tuple val(sampleID), path("*_fastqc.{zip,html}"), emit: quality_stats


    script:
    if (params.workflow == "chipseq" && params.read_type == 'SE')
    """
        [ ! -f  ${sampleID}.fastq.gz ] && ln -s ${fq_reads} ${sampleID}.fastq.gz
        
        fastqc --quiet -t ${task.cpus} ${sampleID}.fastq.gz
    """
    else if (params.workflow == "chipseq" && params.read_type == 'PE')
    """
        [ ! -f  ${sampleID}_1.fastq.gz ] && ln -s ${fq_reads[0]} ${sampleID}_1.fastq.gz
        [ ! -f  ${sampleID}_2.fastq.gz ] && ln -s ${fq_reads[1]} ${sampleID}_2.fastq.gz

        fastqc --quiet -t ${task.cpus} ${sampleID}_1.fastq.gz
        fastqc --quiet -t ${task.cpus} ${sampleID}_2.fastq.gz
    """
    else
    """
        fastqc --quiet -t ${task.cpus} ${fq_reads}
    """

}
