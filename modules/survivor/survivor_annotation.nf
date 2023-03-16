process SURVIVOR_ANNOTATION {
    
    tag "$sampleID"

    cpus 1
    memory 100.GB
    time "00:30:00"

    container 'rocker/tidyverse:4.2.1'    

    publishDir "${params.outdir}", mode:'copy'

    input:
        tuple val(sampleID), file("${sampleID}.ins.bed"),file("${sampleID}.del.bed"),file("${sampleID}.dup.bed"),file("${sampleID}.inv.bed"),file("${sampleID}.tra.bed"), file("${sampleID}.ins.s.bed"), file("${sampleID}.ins.e.bed"), file("${sampleID}.del.s.bed"), file("${sampleID}.del.e.bed"), file("${sampleID}.inv.e.bed"), file("${sampleID}.tra.e.bed"), file("${sampleID}.dup.e.bed"), file("${sampleID}.ins.genes.bed"), file("${sampleID}.del.genes.bed"), file("${sampleID}.inv.genes.bed"), file("${sampleID}.dup.genes.bed"), file("${sampleID}.tra.genes.bed"), file("${sampleID}.ins.exons.bed"), file("${sampleID}.del.exons.bed"), file("${sampleID}.inv.exons.bed"), file("${sampleID}.dup.exons.bed"), file("${sampleID}.tra.exons.bed"), file("${sampleID}.survivor_summary.csv"), file("${sampleID}.merged.overlap.annotated.txt")
    output:
        tuple val(sampleID), file("${sampleID}.survivor_joined_results.csv"), emit: csv
    script:
        """
        /usr/bin/env Rscript ${projectDir}/bin/post_filt.R ${sampleID}
        """
}