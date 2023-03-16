process SURVIVOR_TO_BED {
    tag "$sampleID"

    cpus 1
    memory 100.GB
    time "00:30:00"

    container 'rocker/tidyverse:4.2.1'

    input:
        tuple val(sampleID), file(annot), file(summary)
    output:
        tuple val(sampleID), file("${sampleID}.ins.bed"), file("${sampleID}.del.bed"), file("${sampleID}.dup.bed"), file("${sampleID}.inv.bed"), file("${sampleID}.tra.bed"), emit: sv_beds
    script:
        """
        /usr/bin/env Rscript ${projectDir}/bin/surv_annot_process.R ${sampleID}
        """
}