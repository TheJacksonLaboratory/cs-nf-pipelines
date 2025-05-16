process GATK_REALIGNERTARGETCREATOR {
    // Depricated in GATK4, not recommended. Leaving for historic precedence.
    tag "$sampleID"

    cpus = 12
    memory = 35.GB
    time = '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk3:3.6-0'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.intervals", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*.intervals"), emit: intervals

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /usr/GenomeAnalysisTK.jar \
    -I ${bam} \
    -R ${params.ref_fa} \
    -T RealignerTargetCreator \
    -o ${sampleID}.aligner.intervals \
    -nt $task.cpus \
    --disable_auto_index_creation_and_locking_when_reading_rods
    """
}