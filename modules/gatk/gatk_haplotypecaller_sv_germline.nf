process GATK_HAPLOTYPECALLER_SV_GERMLINE {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time = '10:00:00'

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.*vcf", mode:'copy'

    input:
    tuple val(sampleID), val(meta), file(normal_bam), file(normal_bai)

    output:
    tuple val(sampleID), file("*.*vcf"), emit: vcf
    tuple val(sampleID), file("*.idx"), emit: idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    gatk --java-options "-Xmx${my_mem}G" HaplotypeCaller  \
    -R ${params.ref_fa} \
    -I ${normal_bam} \
    -O ${sampleID}_variants_raw.gvcf \
    -L ${params.target_gatk} \
    -stand-call-conf ${params.call_val} \
    -G StandardAnnotation \
    -G StandardHCAnnotation \
    -G AS_StandardAnnotation \
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
    -ERC GVCF
    """
}