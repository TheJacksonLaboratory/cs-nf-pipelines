process GATK_GETSAMPLENAME {
    tag "$sampleID"

    cpus = 1
    memory = 1.GB
    time = '00:05:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), file(bam), file(bai)

    output:
    tuple val(sampleID), stdout, emit: sample_name

    script:
    """
    mkdir -p tmp
    gatk --java-options "-Djava.io.tmpdir=`pwd`/tmp" GetSampleName \
    -I ${bam} \
    -O sample_name.txt

    cat sample_name.txt
    """
}
