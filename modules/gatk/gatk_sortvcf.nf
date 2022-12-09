process GATK_SORTVCF {

    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '05:30:00'

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'lancet' }", pattern:"*_merged_lancet.vcf.gz", mode:'copy'

    input:
    tuple val(sampleID), path(list)
    val(gvcf)

    output:
    tuple val(sampleID), file("*.vcf"), file("*.idx"), emit: vcf_idx, optional: true
    tuple val(sampleID), file("*_merged_lancet.vcf.gz"), emit: lancet_vcf, optional: true

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    inputs = list.collect { "-I $it" }.join(' ')

    lancet_check = (!!(list[0] =~ /lancet/))

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

    if [ $lancet_check = true ]; then
        mv ${sampleID}_merged.${output_suffix} ${sampleID}_merged_lancet.vcf
        bgzip ${sampleID}_merged_lancet.vcf
    fi
    """
}


    //  ; then

    // fi