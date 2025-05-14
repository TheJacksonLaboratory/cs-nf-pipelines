process CLUMPIFY {
    tag "$sampleID"

    cpus 6
    memory 200.GB
    time '48:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bbmap:39.06--h92535d8_0'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.txt", mode:'copy'

    input:
    tuple val(sampleID), path(fq_reads)

    output:
    tuple val(sampleID), path("${sampleID}.clumpy.R*.fastq.gz"), emit: clumpy_fastq
    tuple val(sampleID), path("*log.txt"), emit: clumpy_log

    script:
    if (params.read_type == "SE")
    """
    testformat.sh ${fq_reads[0]} > fastq_format.txt

    if grep -q 'illumina' fastq_format.txt ; then qual=64; else qual=33; fi

    clumpify.sh in=${fq_reads[0]} out=${sampleID}.clumpy.R1.fastq.gz tmpdir=./ usetmpdir=t dedupe=t qin=\${qual} -Xmx199g &> ${sampleID}_clumpy_log.txt
    """
    else
    """
    testformat.sh ${fq_reads[0]} > fastq_format.txt
    
    if grep -q 'illumina' fastq_format.txt ; then qual=64; else qual=33; fi

    clumpify.sh in=${fq_reads[0]} in2=${fq_reads[1]} out=${sampleID}.clumpy.R1.fastq.gz out2=${sampleID}.clumpy.R2.fastq.gz tmpdir=./ usetmpdir=t dedupe=t qin=\${qual} -Xmx199g &> ${sampleID}_clumpy_log.txt
    """
}
