process SURVIVOR_INEXON {
    tag "$sampleID"

    cpus 1
    memory 20.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/pysam:0.15.2--py36h02877da_7'

    publishDir "${params.pubdir}/${sampleID}", mode:'copy', enabled: params.data_type == 'ont' ? false : true

    input:
        tuple val(sampleID), file(survivor_vcf), file("ins.exons.bed"), file("del.exons.bed"), file("dup.exons.bed"), file("tra.exons.bed"), file("inv.exons.bed")
    output:
        tuple val(sampleID), file("${sampleID}_*_struct_var.vcf"), emit: vcf
    script:

    if (params.data_type == "pacbio")
        """
        /usr/bin/env python ${projectDir}/bin/germline_sv/annot_vcf_with_exon.py -v ${survivor_vcf} \
            -i ins.exons.bed -d del.exons.bed \
            -u dup.exons.bed -t tra.exons.bed -n inv.exons.bed \
            -o ${sampleID}_PACBIO_PS_struct_var.vcf
        """
    else if (params.data_type == "illumina")
        """
        /usr/bin/env python ${projectDir}/bin/germline_sv/annot_vcf_with_exon.py -v ${survivor_vcf} \
            -i ins.exons.bed -d del.exons.bed \
            -u dup.exons.bed -t tra.exons.bed -n inv.exons.bed \
            -o ${sampleID}_ILLUMINA_DLM_struct_var.vcf
        """
    else if (params.data_type == "ont")
        """
        /usr/bin/env python ${projectDir}/bin/germline_sv/annot_vcf_with_exon.py -v ${survivor_vcf} \
            -i ins.exons.bed -d del.exons.bed \
            -u dup.exons.bed -t tra.exons.bed -n inv.exons.bed \
            -o ${sampleID}_PACBIO_NS_struct_var.vcf
        """           
}
