process GRIDSS_CALLING {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time = '10:00:00'

    container 'quay.io/jaxcompsci/gridss:2.13.2-2_ln'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gridss' }", pattern: "*_gridss_sv.vcf.gz", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), val(meta), file(normal_bam), file(normal_bai), val(normal_name), file(tumor_bam), file(tumor_bai), val(tumor_name), val(gridss_assembled)

    output:
    tuple val(sampleID), path('*_gridss_sv.vcf.gz'), val(meta), val(normal_name), val(tumor_name), emit: gridss_vcf

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]+'g'
    
    output_dir = 'gridss_call/'

    """
    # https://github.com/umccr/gridss-purple-linx-nf
    # Create shadow directory with file symlinks of GRIDSS 'workingdir' to prevent NF cache invalidation (resume related)
    # NOTE: for reasons that elude me, NF doesn't always stage in the workingdir; remove if it is present
    mkdir -p "${output_dir}"
    lndir \$(readlink -f "${gridss_assembled}/") "${output_dir}/"
    if [[ -L "${gridss_assembled.name}" ]]; then
        rm "${gridss_assembled}"
    fi

    gridss \
    --jvmheap "${my_mem}" \
    --steps call \
    --reference "${params.ref_fa_indices}" \
    --jar /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    --threads ${task.cpus} \
    --workingdir "${output_dir}/work/" \
    --assembly "${output_dir}/${sampleID}.gridssassembly.bam" \
    --output "${sampleID}_gridss_sv.vcf.gz" \
    --picardoptions VALIDATION_STRINGENCY=LENIENT \
    ${normal_bam} ${tumor_bam}
    """
}