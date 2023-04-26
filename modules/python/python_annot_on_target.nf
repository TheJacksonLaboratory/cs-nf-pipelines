process PYTHON_ANNOT_ON_TARGET {

    tag "$sampleID"

    cpus 1
    memory 20.GB
    time "00:30:00"

    container 'quay.io/biocontainers/pysam:0.15.2--py36h02877da_7'

    input:
        tuple val(sampleID), path(vcf)
    output:
        tuple val(sampleID), path("${sampleID}_ONT_NS_struct_var.vcf"), emit: vcf
    
    publishDir "${params.pubdir}", mode:'copy'
    
    script:

    if (params.workflow == "ont" && params.targ_chr && params.targ_start && params.targ_end)
        """
        /usr/bin/env python ${projectDir}/bin/annot_vcf_with_on_target.py \
            -v ${vcf} \
            -c ${params.targ_chr} \
            -s ${params.targ_start} \
            -e ${params.targ_end} \
            -o ${sampleID}_ONT_NS_struct_var.vcf
        """
    else if (params.workflow == "ont")
        """
        /usr/bin/env python ${projectDir}/bin/annot_vcf_with_on_target.py \
            -v ${vcf} \
            -c NA \
            -s NA \
            -e NA \
            -o ${sampleID}_ONT_NS_struct_var.vcf
        """
    else error "module relies on script that currently only supports ONT"
}