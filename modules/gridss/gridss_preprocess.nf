process GRIDSS_PREPROCESS {
    tag "$meta.patient"

    cpus = 4
    memory = 15.GB
    time = '10:00:00'

    container 'quay.io/jaxcompsci/gridss:2.13.2-2_ln'

    input:
    tuple val(sampleID), val(meta), file(normal_bam), file(normal_bai), val(normal_name), file(tumor_bam), file(tumor_bai), val(tumor_name)

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