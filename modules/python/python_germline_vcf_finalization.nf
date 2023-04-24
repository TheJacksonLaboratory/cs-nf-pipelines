process GERMLINE_VCF_FINALIZATION {
    tag "$sampleID"

    cpus 1
    memory { 5.GB * task.attempt }
    time {1.hour * task.attempt}
    errorStrategy 'retry'
    maxRetries 1

    container 'quay.io/jaxcompsci/py3_perl_pylibs:v2'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'vcf' }", pattern: "*final.vcf", mode:'copy'

    input:
    tuple val(sampleID), file(vcf)
    val(filtered)

    output:
    tuple val(sampleID), file("*final.vcf"), emit: vcf

    script:

    output_suffix = filtered == 'filtered' ?  'filtered' : 'unfiltered'

    """
    python \
    ${projectDir}/bin/pta/annotate_id.py \
    ${vcf} \
    ${sampleID}_germline_vep_cosmic_cancerResitMut_annotated_id.vcf

    python \
    ${projectDir}/bin/pta/rename_csq_vcf.py \
    ${sampleID}_germline_vep_cosmic_cancerResitMut_annotated_id.vcf \
    ${sampleID}_germline_snv_indel_annotated_${output_suffix}_final.vcf

    """
}