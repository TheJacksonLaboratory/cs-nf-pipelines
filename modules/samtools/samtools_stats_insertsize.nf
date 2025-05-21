process SAMTOOLS_STATS_INSERTSIZE {
    tag "$sampleID"

    cpus 8
    memory 1.GB
    time '02:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*insert_size.txt", mode:'copy'

    input:
        tuple val(sampleID), val(meta), path(bam), path(bai), val(read_ID)

    output:
        tuple val(sampleID), env(read_length), env(insert_size), emit: read_length_insert_size
        file("*insert_size.txt")

    script:
    """
    samtools stats --insert-size 8000  ${bam} --threads ${task.cpus} | grep ^SN | cut -f 2- > ${sampleID}_insert_size.txt
    read_length=`grep "maximum length" ${sampleID}_insert_size.txt | cut -d ':' -f2 | tr -d " \\t\\n\\r"`
    insert_size=`grep "insert size average" ${sampleID}_insert_size.txt | cut -d ':' -f2 | tr -d " \\t\\n\\r"`
    """
}
