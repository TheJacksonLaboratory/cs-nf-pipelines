process SURVIVOR_BED_INTERSECT {
    tag "$sampleID"

    cpus 1
    memory 100.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/bedtools-sv_refs:2.30.0--refv0.2.0'

    input:
        tuple val(sampleID), file("${sampleID}.ins.bed"),file("${sampleID}.del.bed"),file("${sampleID}.dup.bed"),file("${sampleID}.inv.bed"),file("${sampleID}.tra.bed")
    output:
        tuple val(sampleID), file("${sampleID}.ins.c.bed"), file("${sampleID}.del.c.bed"), file("${sampleID}.inv.c.bed"), file("${sampleID}.ins.genes.bed"), file("${sampleID}.del.genes.bed"), file("${sampleID}.inv.genes.bed"), file("${sampleID}.dup.genes.bed"), file("${sampleID}.tra.genes.bed"), file("${sampleID}.ins.exons.bed"), file("${sampleID}.del.exons.bed"), file("${sampleID}.inv.exons.bed"), file("${sampleID}.dup.exons.bed"), file("${sampleID}.tra.exons.bed"), file("${sampleID}.ins.regulatory.bed"), file("${sampleID}.del.regulatory.bed"), file("${sampleID}.inv.regulatory.bed"), file("${sampleID}.dup.regulatory.bed"), file("${sampleID}.tra.regulatory.bed"), emit: intersected_beds
        tuple val(sampleID), file("${sampleID}.ins.exons.bed"), file("${sampleID}.del.exons.bed"), file("${sampleID}.inv.exons.bed"), file("${sampleID}.dup.exons.bed"), file("${sampleID}.tra.exons.bed"), emit: intersected_exons
    script:
        """
        bedtools intersect -a ${sampleID}.ins.bed -b ${params.sv_ins_ref} -f 0.2 > ${sampleID}.ins.c.bed
        bedtools intersect -a ${sampleID}.del.bed -b ${params.sv_del_ref} -f 0.2 > ${sampleID}.del.c.bed
        bedtools intersect -a ${sampleID}.inv.bed -b ${params.sv_inv_ref} -f 0.2 > ${sampleID}.inv.c.bed

        bedtools intersect -a ${sampleID}.ins.bed -b ${params.genes_bed} -wa -wb > ${sampleID}.ins.genes.bed
        bedtools intersect -a ${sampleID}.del.bed -b ${params.genes_bed} -wa -wb > ${sampleID}.del.genes.bed
        bedtools intersect -a ${sampleID}.inv.bed -b ${params.genes_bed} -wa -wb > ${sampleID}.inv.genes.bed
        bedtools intersect -a ${sampleID}.dup.bed -b ${params.genes_bed} -wa -wb > ${sampleID}.dup.genes.bed
        bedtools intersect -a ${sampleID}.tra.bed -b ${params.genes_bed} -wa -wb > ${sampleID}.tra.genes.bed

        bedtools intersect -a ${sampleID}.ins.bed -b ${params.exons_bed} -wa -wb > ${sampleID}.ins.exons.bed
        bedtools intersect -a  ${sampleID}.del.bed -b ${params.exons_bed} -wa -wb > ${sampleID}.del.exons.bed
        bedtools intersect -a  ${sampleID}.inv.bed -b ${params.exons_bed} -wa -wb > ${sampleID}.inv.exons.bed
        bedtools intersect -a  ${sampleID}.dup.bed -b ${params.exons_bed} -wa -wb > ${sampleID}.dup.exons.bed
        bedtools intersect -a  ${sampleID}.tra.bed -b ${params.exons_bed} -wa -wb > ${sampleID}.tra.exons.bed

        bedtools intersect -a ${sampleID}.ins.bed -b ${params.reg_ref} -wa -wb > ${sampleID}.ins.regulatory.bed
        bedtools intersect -a ${sampleID}.del.bed -b ${params.reg_ref} -wa -wb > ${sampleID}.del.regulatory.bed
        bedtools intersect -a ${sampleID}.inv.bed -b ${params.reg_ref} -wa -wb > ${sampleID}.inv.regulatory.bed
        bedtools intersect -a ${sampleID}.dup.bed -b ${params.reg_ref} -wa -wb > ${sampleID}.dup.regulatory.bed
        bedtools intersect -a ${sampleID}.tra.bed -b ${params.reg_ref} -wa -wb > ${sampleID}.tra.regulatory.bed
        """
}
