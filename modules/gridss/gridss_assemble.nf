process GRIDSS_ASSEMBLE {
    tag "$sampleID"

    cpus = 4
    memory = 100.GB
    time = '20:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/gridss:2.13.2-3'

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name), val(gridss_preprocessed)

    output:
    tuple val(sampleID), path('gridss_assemble/'), emit: gridss_assembly

    script:
    String my_mem = (task.memory-5.GB).toString()
    String my_other_mem = (task.memory-30.GB).toString()
    heap_mem =  my_mem[0..-4]+'g'
    other_mem = my_other_mem[0..-4]+'g'

    output_dir = 'gridss_assemble/'

    """
    # https://github.com/umccr/gridss-purple-linx-nf
    # Create shadow directory with file symlinks of GRIDSS 'workingdir' to prevent NF cache invalidation (resume related)
    # NOTE: for reasons that elude me, NF doesn't always stage in the workingdir; remove if it is present
    mkdir -p "${output_dir}/work/"
    lndir \$(readlink -f "${gridss_preprocessed}/") "${output_dir}/work"
    if [[ -L "${gridss_preprocessed.name}" ]]; then
        rm "${gridss_preprocessed}"
    fi

    gridss \
    --jvmheap "${heap_mem}" \
    --otherjvmheap "${other_mem}" \
    --steps assemble \
    --reference "${params.combined_reference_set}" \
    --jar /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    --threads ${task.cpus} \
    --workingdir "${output_dir}/work/" \
    --assembly ${output_dir}/${sampleID}.gridssassembly.bam \
    --picardoptions VALIDATION_STRINGENCY=LENIENT \
    ${normal_bam} ${tumor_bam}
    """
}
