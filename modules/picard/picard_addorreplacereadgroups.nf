process PICARD_ADDORREPLACEREADGROUPS {
    tag "$sampleID"

    cpus 1
    memory 45.GB
    time '16:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*.bam", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(read_groups), file(bam)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam
    tuple val(sampleID), file("*.bai"), emit: bai

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp AddOrReplaceReadGroups \
    INPUT=${bam} \
    OUTPUT=${sampleID}_genome_bam_with_read_groups.bam \
    SORT_ORDER=coordinate \
    TMP_DIR=`pwd`/tmp \
    \$(cat $read_groups) \
    CREATE_INDEX=true
    """
}
