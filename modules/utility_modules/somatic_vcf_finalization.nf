process SOMATIC_VCF_FINALIZATION {
    tag "$sampleID"

    cpus 1
    memory { 5.GB * task.attempt }
    time {1.hour * task.attempt}
    errorStrategy 'retry'
    maxRetries 1

    container 'quay.io/jaxcompsci/py3_perl_pylibs:v2'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'vcf' }", pattern: "*final.*", mode:'copy'

    input:
    tuple val(sampleID), file(vcf), val(meta), val(normal_name), val(tumor_name)
    val(filtered)

    output:
    tuple val(sampleID), file("*final.vcf"), emit: vcf

    script:

    output_suffix = filtered == 'filtered' ?  'filtered' : 'unfiltered'

    """
    python \
    ${projectDir}/bin/sv/annotate_id.py \
    ${vcf} \
    ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_id.vcf

    python \
    ${projectDir}/bin/sv/rename_csq_vcf.py \
    ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_id.vcf \
    ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_${output_suffix}_csq.vcf

    python \
    ${projectDir}/bin/sv/make_main_vcf.py \
    ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_${output_suffix}_csq.vcf \
    ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_${output_suffix}_final.vcf

    python \
    ${projectDir}/bin/sv/make_txt.py \
    --vcf ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_${output_suffix}_final.vcf \
    --txt ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_${output_suffix}_final.txt \
    --tumor ${tumor_name} \
    --normal ${normal_name} 

    python \
    ${projectDir}/bin/sv/make_maf.py \
    --vcf ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_${output_suffix}_final.vcf \
    --maf ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_${output_suffix}_final.maf \
    --library WGS \
    --vep-version GRCh38 \
    --tumor ${tumor_name} \
    --normal ${normal_name} \
    --ensembl-entrez ${params.ensembl_entrez}

    """
}