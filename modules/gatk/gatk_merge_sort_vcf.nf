process GATK_SORTMERGEVCF {

    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '05:30:00'

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), path(list)

    output:
    tuple val(sampleID), file("*.vcf"), file("*.idx"), emit: vcf_idx
    tuple val(sampleID), file("*merged_sort.vcf"), emit: merge_sort_vcf, 

    script:

    inputs = list.collect { "-I $it" }.join(' ')

    """
    gatk --java-options "-Xmx${my_mem}G" SortVcf  \
        -SD ${params.ref_fa_dict} \
        ${inputs} \
        -O ${sampleID}_merged_sort.vcf
    """
}