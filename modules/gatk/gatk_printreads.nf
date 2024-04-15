process GATK_PRINTREADS {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.4.0.0'

    input:
    tuple val(sampleID), val(meta), path(bam), path(bai)

    output:
    tuple val(sampleID), val(meta), file("*.bam"), file("*.bai"), emit: bam_bai

    script:
    //Estimate somatic variants using Mutect2
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir tmp
    gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=`pwd`/tmp" PrintReads \
    -L ${params.target_gatk} \
    -I ${bam} \
    -O ${bam.baseName}.targetOnly.bam
    """
}
