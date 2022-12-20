
process GRIDSS_SOMATIC_FILTER {
    tag "$sampleID"

    cpus = 1
    memory = 5.GB
    time = '01:00:00'
    errorStrategy 'ignore'
    
    container 'quay.io/jaxcompsci/gridss:2.13.2-2_ln'

    stageInMode = 'copy'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gridss' }", pattern: "*_gridss_sv_somaticFiltered.vcf", mode:'copy'

    input:
    tuple val(sampleID), path(vcf)
    path(gridss_pon)

    output:
    tuple val(sampleID), path('*_gridss_sv_somaticFiltered.vcf'), emit: gridss_filtered_vcf

    script:

    """
    tar -zxvf ${gridss_pon}

    Rscript /opt/gridss/gridss_somatic_filter \
    --ref ${params.gridss_somatic_filter_ref} \
    --input ${vcf} \
    --output ${sampleID}_gridss_sv_somaticFiltered.vcf \
    --tumourordinal 2 \
    --plotdir . \
    --scriptdir /opt/gridss/ \
    --configdir /opt/gridss/ \
    --pondir pon/
    """
        
    stub:
    """
    touch ${sampleID}_gridss_sv_somaticFiltered.vcf
    """
}