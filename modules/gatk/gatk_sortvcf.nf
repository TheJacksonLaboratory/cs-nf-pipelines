process GATK_SORTVCF {

    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '05:30:00'

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), path(list)
    val(gvcf)

    output:
    tuple val(sampleID), file("*.vcf"), file("*.idx"), emit: vcf_idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    inputs = list.collect { "-I $it" }.join(' ')

    if (gvcf=='gvcf'){
        output_suffix='g.vcf'
    }
    else{
        output_suffix='vcf'
    }

    """
    gatk --java-options "-Xmx${my_mem}G" SortVcf  \
        -SD ${params.ref_fa_dict} \
        ${inputs} \
        -O ${sampleID}_merged.${output_suffix}
    """
}