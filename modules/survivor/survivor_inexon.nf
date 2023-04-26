process SURVIVOR_INEXON {

    tag "$sampleID"

    cpus 1
    memory 20.GB
    time "00:30:00"

    container 'quay.io/biocontainers/pysam:0.15.2--py36h02877da_7'

    publishDir "${params.pubdir}", mode:'copy', enabled: params.workflow == 'ont' ? false : true

    input:
        tuple val(sampleID), file(survivor_vcf), file("ins.exons.bed"), file("del.exons.bed"), file("dup.exons.bed"), file("tra.exons.bed"), file("inv.exons.bed")
    output:
        tuple val(sampleID), file("${sampleID}_*_struct_var.vcf"), emit: vcf
    script:

    if (params.workflow == "pacbio")
        """
        /usr/bin/env python ${projectDir}/bin/annot_vcf_with_exon.py -v ${survivor_vcf} \
            -i ins.exons.bed -d del.exons.bed \
            -u dup.exons.bed -t tra.exons.bed -n inv.exons.bed \
            -o ${sampleID}_PACBIO_PS_struct_var.vcf
        """
    else if (params.workflow == "illumina")
        """
        /usr/bin/env python ${projectDir}/bin/annot_vcf_with_exon.py -v ${survivor_vcf} \
            -i ins.exons.bed -d del.exons.bed \
            -u dup.exons.bed -t tra.exons.bed -n inv.exons.bed \
            -o ${sampleID}_ILLUMINA_BDLM_struct_var.vcf
        """
    else if (params.workflow == "ont")
        """
        /usr/bin/env python ${projectDir}/bin/annot_vcf_with_exon.py -v ${survivor_vcf} \
            -i ins.exons.bed -d del.exons.bed \
            -u dup.exons.bed -t tra.exons.bed -n inv.exons.bed \
            -o ${sampleID}_PACBIO_NS_struct_var.vcf
        """           
}