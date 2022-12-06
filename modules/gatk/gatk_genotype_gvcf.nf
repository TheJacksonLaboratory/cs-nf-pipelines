process GATK_GENOTYPE_GVCF {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '01:30:00'

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.*vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(vcf), file(vcf_index), path(interval), val(index)

    output:
    tuple val(sampleID), file("*.*vcf"), file("*.idx"), path(interval), val(index), emit: vcf_idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    gatk --java-options "-Xmx${my_mem}G" GenotypeGVCFs  \
    -R ${params.ref_fa} \
    -V ${vcf} \
    -O ${sampleID}_${index}_genotypedGVCFs.vcf \
    -L ${interval}
    """
}