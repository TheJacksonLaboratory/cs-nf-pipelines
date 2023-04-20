process SOMATIC_VCF_FINALIZATION {
    tag "$sampleID"

    cpus 1
    memory 50.GB
    time 1.hour
    
    container 'quay.io/jaxcompsci/py3_perl_pylibs:v2'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'vcf' }", pattern: "*final.*", mode:'copy'
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'vcf' }", pattern: "*supplemental.vcf", mode:'copy'

    input:
    tuple val(sampleID), file(vcf), val(meta), val(normal_name), val(tumor_name)
    val(filtered)

    output:
    tuple val(sampleID), file("*final.vcf"), emit: vcf
    tuple val(sampleID), file("*final.txt"), emit: txt
    tuple val(sampleID), file("*final.maf"), emit: maf
    tuple val(sampleID), file("*supplemental.vcf"), emit: supp_vcf

    script:

    output_suffix = filtered == 'filtered' ?  'filtered' : 'unfiltered'

    """
    python \
    ${projectDir}/bin/pta/annotate_id.py \
    ${vcf} \
    ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_id.vcf

    python \
    ${projectDir}/bin/pta/rename_csq_vcf.py \
    ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_id.vcf \
    ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_${output_suffix}_supplemental.vcf

    python \
    ${projectDir}/bin/pta/make_main_vcf.py \
    ${sampleID}_somatic_vep_cosmic_cancerResitMut_annotated_${output_suffix}_supplemental.vcf \
    ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_final.vcf

    python \
    ${projectDir}/bin/pta/make_txt.py \
    --vcf ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_final.vcf \
    --txt ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_final.txt \
    --tumor ${tumor_name} \
    --normal ${normal_name} 

    python \
    ${projectDir}/bin/pta/make_maf.py \
    --vcf ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_final.vcf \
    --maf ${sampleID}_somatic_snv_indel_annotated_${output_suffix}_final.maf \
    --library WGS \
    --vep-version GRCh38 \
    --tumor ${tumor_name} \
    --normal ${normal_name} \
    --ensembl-entrez ${params.ensembl_entrez}

    """
}