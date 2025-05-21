process BCFTOOLS_QUERY_ASCAT {
    tag "$sampleID"
    
    cpus 1
    memory 8.GB
    time '01:00:00'
    errorStrategy 'finish'

    container 'quay.io/jaxcompsci/gtc2vcf_with_tools:v2'
    publishDir "${params.pubdir}/${sampleID}", mode: 'copy'

    input:
    tuple val(sampleID), val(meta), path(bcf), path(csi), path(vcf), path(tsv)

    output:
    tuple val(sampleID), val(meta), path('*_convert.BAF'), path('*_convert.LRR'), emit: baf_lrr

    script:
    """
    (bcftools query -l ${sampleID}_convert.bcf | awk 'BEGIN{printf("\\tCHROM\\tPOS");} {printf("\\t%s",\$1);} END{printf("\\n");}'  &&  bcftools query -f '%ID\\t%CHROM\\t%POS[\\t%BAF]\\n' ${sampleID}_convert.bcf) > ${sampleID}_convert.BAF
    
    (bcftools query -l ${sampleID}_convert.bcf | awk 'BEGIN{printf("\\tCHROM\\tPOS");} {printf("\\t%s",\$1);} END{printf("\\n");}'  &&  bcftools query -f '%ID\\t%CHROM\\t%POS[\\t%LRR]\\n' ${sampleID}_convert.bcf) > ${sampleID}_convert.LRR
    
    sed -i s/chr// ${sampleID}_convert.BAF
    sed -i s/chr// ${sampleID}_convert.LRR
    """
}
