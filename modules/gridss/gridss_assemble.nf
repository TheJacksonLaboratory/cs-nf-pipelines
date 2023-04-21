process GRIDSS_ASSEMBLE {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time = '10:00:00'

    container 'quay.io/jaxcompsci/gridss:2.13.2-2_ln'

    input:
    tuple val(sampleID), val(meta), file(normal_bam), file(normal_bai), val(normal_name), file(tumor_bam), file(tumor_bai), val(tumor_name), val(gridss_preprocessed)

    output:
    tuple val(sampleID), path('gridss_assemble/'), emit: gridss_assembly

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]+'g'
    
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
    --jvmheap "${my_mem}" \
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