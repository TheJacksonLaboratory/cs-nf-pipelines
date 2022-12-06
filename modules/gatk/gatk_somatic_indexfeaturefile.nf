process GATK_INDEX_FEATURE_FILE_SV_SOMATIC {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '10:00:00'

    container 'broadinstitute/gatk:4.2.4.1'
    input:
    tuple val(sampleID), val(meta), file(vcf_compressed)

    output:
    tuple val(sampleID), file("*.*vcf.gz"), emit: vcf_index

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    gatk --java-options "-Xmx${my_mem}G" IndexFeatureFile  \
    --input ${sampleID}.vcf.gz 
    """
}