process GATK_INDELREALIGNER{
    // Deprecated in GATK4, not recommended. Leaving for historic precedence.
    tag "$sampleID"

    cpus = 1
    memory = 35.GB
    time = '08:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    // Command Depricated in GATK 4
    container 'broadinstitute/gatk3:3.6-0'

    // save if mouse and wgs or save if keep intermediate
    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*.bam", mode:'copy', enabled: params.gen_org=='mouse' && params.workflow=='wgs' ? true : params.keep_intermediate

    input:
    tuple val(sampleID), file(bam)
    tuple val(sampleID), file(intervals)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam
    tuple val(sampleID), file("*.bai"), emit: bai

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /usr/GenomeAnalysisTK.jar \
    -I ${bam} \
    -R ${params.ref_fa} \
    -T IndelRealigner \
    -targetIntervals ${intervals} \
    -o ${sampleID}_realigned.bam \
    --disable_auto_index_creation_and_locking_when_reading_rods
    """
}