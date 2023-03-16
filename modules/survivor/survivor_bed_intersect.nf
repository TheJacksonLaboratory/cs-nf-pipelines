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
        /usr/bin/env bash ${projectDir}/bin/intersect_beds.sh ${sampleID}
        """
}