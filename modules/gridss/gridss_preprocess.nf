process GRIDSS_PREPROCESS {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time = '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gridss:2.13.2-3'

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name)

    output:
    tuple val(sampleID), path('gridss_preprocess/'), emit: gridss_preproc

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]+'g'
    """
    # https://github.com/umccr/gridss-purple-linx-nf
    gridss \
    --jvmheap "${my_mem}" \
    --steps preprocess \
    --reference "${params.combined_reference_set}" \
    --jar /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    --threads ${task.cpus} \
    --workingdir gridss_preprocess/ \
    --picardoptions VALIDATION_STRINGENCY=LENIENT \
    ${normal_bam} ${tumor_bam}
    """
}
