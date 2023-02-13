process GRIDSS_CHROM_FILTER {
    tag "$sampleID"

    cpus = 1
    memory = 1.GB
    time = '01:00:00'
    
    container 'quay.io/jaxcompsci/internal_tools:v1.0'

    stageInMode = 'copy'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gridss' }", pattern: "*_gridss_sv_unfiltered_chroms.vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(vcf), val(meta), val(normal_name), val(tumor_name)
    val(chroms)

    output:
    tuple val(sampleID), path('*_gridss_sv_unfiltered_chroms.vcf'), val(meta), val(normal_name), val(tumor_name), emit: gridss_chrom_vcf
    
    script:
    chrom_list = chroms.collect { "$it" }.join(' ')

    """
    python ${projectDir}/bin/sv/filter_vcf.py \
    --vcf-file ${vcf} \
    --output ${sampleID}_gridss_sv_unfiltered_chroms.vcf \
    --chroms ${chrom_list}
    """
}
