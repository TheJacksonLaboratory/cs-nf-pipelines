process SURVIVOR_ANNOTATION {
    tag "$sampleID"

    cpus 1
    memory 100.GB
    time "01:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'rocker/tidyverse:4.2.1'    

    publishDir "${params.pubdir}/${sampleID}", mode:'copy', enabled: params.data_type == 'ont' ? false : true

    input:
        tuple val(sampleID), file("${sampleID}.ins.bed"), file("${sampleID}.del.bed"), file("${sampleID}.dup.bed"), file("${sampleID}.inv.bed"), file("${sampleID}.tra.bed"), file("${sampleID}.ins.c.bed"), file("${sampleID}.del.c.bed"), file("${sampleID}.inv.c.bed"), file("${sampleID}.ins.genes.bed"), file("${sampleID}.del.genes.bed"), file("${sampleID}.inv.genes.bed"), file("${sampleID}.dup.genes.bed"), file("${sampleID}.tra.genes.bed"), file("${sampleID}.ins.exons.bed"), file("${sampleID}.del.exons.bed"), file("${sampleID}.inv.exons.bed"), file("${sampleID}.dup.exons.bed"), file("${sampleID}.tra.exons.bed"), file("${sampleID}.ins.regulatory.bed"), file("${sampleID}.del.regulatory.bed"), file("${sampleID}.inv.regulatory.bed"), file("${sampleID}.dup.regulatory.bed"), file("${sampleID}.tra.regulatory.bed"), file("${sampleID}.survivor_summary.csv"), file("${sampleID}.merged.overlap.annotated.txt")
    output:
        tuple val(sampleID), file("*survivor_joined_results.csv"), emit: csv
    script:
        """
        /usr/bin/env Rscript ${projectDir}/bin/germline_sv/summarize_intersections.R ${sampleID}
        """
}
