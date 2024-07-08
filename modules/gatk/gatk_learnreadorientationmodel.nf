process GATK_LEARNREADORIENTATIONMODEL {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.4.0.0'

    input:
    tuple val(sampleID), path(f1r2)

    output:
    tuple val(sampleID), file("*readOrientationModel.tar.gz"), emit: model_file

    script:
    //Estimate somatic variants using Mutect2
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" LearnReadOrientationModel \
    -I ${f1r2} \
    -O ${sampleID}.readOrientationModel.tar.gz
    """
}
