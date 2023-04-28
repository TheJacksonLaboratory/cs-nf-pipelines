process SURVIVOR_BED_INTERSECT {
    tag "$sampleID"

    cpus 1
    memory 100.GB
    time "00:30:00"

    container 'quay.io/jaxcompsci/bedtools-sv_refs:2.30.0--hc088bd4_0'

    input:
        tuple val(sampleID), file("${sampleID}.ins.bed"),file("${sampleID}.del.bed"),file("${sampleID}.dup.bed"),file("${sampleID}.inv.bed"),file("${sampleID}.tra.bed")
    output:
        tuple val(sampleID), file("${sampleID}.ins.s.bed"), file("${sampleID}.ins.e.bed"), file("${sampleID}.del.s.bed"), file("${sampleID}.del.e.bed"), file("${sampleID}.inv.e.bed"), file("${sampleID}.tra.e.bed"), file("${sampleID}.dup.e.bed"), file("${sampleID}.ins.genes.bed"), file("${sampleID}.del.genes.bed"), file("${sampleID}.inv.genes.bed"), file("${sampleID}.dup.genes.bed"), file("${sampleID}.tra.genes.bed"), file("${sampleID}.ins.exons.bed"), file("${sampleID}.del.exons.bed"), file("${sampleID}.inv.exons.bed"), file("${sampleID}.dup.exons.bed"), file("${sampleID}.tra.exons.bed"), emit: intersected_beds
        tuple val(sampleID), file("${sampleID}.ins.exons.bed"), file("${sampleID}.del.exons.bed"), file("${sampleID}.inv.exons.bed"), file("${sampleID}.dup.exons.bed"), file("${sampleID}.tra.exons.bed"), emit: intersected_exons
    script:
        """
        bedtools window -w 100 -a ${sampleID}.ins.bed -b ${params.sanger_ins} > ${sampleID}.ins.s.bed
        bedtools window -w 100 -a ${sampleID}.ins.bed -b ${params.dgva_ins} > ${sampleID}.ins.e.bed
        bedtools window -w 100 -a ${sampleID}.del.bed -b ${params.sanger_del} > ${sampleID}.del.s.bed
        bedtools window -w 100 -a ${sampleID}.del.bed -b ${params.dgva_del} > ${sampleID}.del.e.bed
        bedtools window -w 100 -a ${sampleID}.inv.bed -b ${params.dgva_inv} > ${sampleID}.inv.e.bed
        bedtools window -w 100 -a ${sampleID}.dup.bed -b ${params.dgva_dup} > ${sampleID}.dup.e.bed
        bedtools window -w 100 -a ${sampleID}.tra.bed -b ${params.dgva_tra} > ${sampleID}.tra.e.bed

        bedtools intersect -a ${sampleID}.ins.bed -b ${params.genes_bed} -wa -wb > ${sampleID}.ins.genes.bed
        bedtools intersect -a ${sampleID}.del.bed -b ${params.genes_bed} -wa -wb > ${sampleID}.del.genes.bed
        bedtools intersect -a ${sampleID}.inv.bed -b ${params.genes_bed} -wa -wb > ${sampleID}.inv.genes.bed
        bedtools intersect -a ${sampleID}.dup.bed -b ${params.genes_bed} -wa -wb > ${sampleID}.dup.genes.bed
        bedtools intersect -a ${sampleID}.tra.bed -b ${params.genes_bed} -wa -wb > ${sampleID}.tra.genes.bed

        bedtools intersect -a ${sampleID}.ins.bed -b ${params.exons_bed} -wa -wb | \
            cut -f 1,2,3,4,5,6,7,8,10,12,14,16 > ${sampleID}.ins.exons.bed
        bedtools intersect -a  ${sampleID}.del.bed -b ${params.exons_bed} -wa -wb | \
            cut -f 1,2,3,4,5,6,7,8,10,12,14,16 > ${sampleID}.del.exons.bed
        bedtools intersect -a  ${sampleID}.inv.bed -b ${params.exons_bed} -wa -wb | \
            cut -f 1,2,3,4,5,6,7,8,10,12,14,16 > ${sampleID}.inv.exons.bed
        bedtools intersect -a  ${sampleID}.dup.bed -b ${params.exons_bed} -wa -wb | \
            cut -f 1,2,3,4,5,6,7,8,10,12,14,16 > ${sampleID}.dup.exons.bed
        bedtools intersect -a  ${sampleID}.tra.bed -b ${params.exons_bed} -wa -wb | \
            cut -f 1,2,3,4,5,6,7,8,10,12,14,16 > ${sampleID}.tra.exons.bed

        bash ${projectDir}/bin/sed_unquote.sh ${sampleID}.ins.exons.bed
        bash ${projectDir}/bin/sed_unquote.sh ${sampleID}.del.exons.bed
        bash ${projectDir}/bin/sed_unquote.sh ${sampleID}.inv.exons.bed
        bash ${projectDir}/bin/sed_unquote.sh ${sampleID}.dup.exons.bed
        bash ${projectDir}/bin/sed_unquote.sh ${sampleID}.tra.exons.bed
        """
}