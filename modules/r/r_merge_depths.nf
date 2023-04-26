process R_MERGE_DEPTHS {
    
    tag "$sampleID"

    cpus 1
    memory 100.GB
    time "00:30:00"

    container 'rocker/tidyverse:4.2.1'    

    publishDir "${params.pubdir}", mode:'copy', pattern: "${sampleID}.survivor_joined_results_depths.csv"

    input:
        tuple val(sampleID), path(nanosv_depths), path(sniffles_depths), path(survivor_ids), path(summary_table)
    output:
        tuple val(sampleID), file("${sampleID}.survivor_joined_results_depths.csv"), emit: csv
        tuple val(sampleID), file("${sampleID}.merged_depths.bed"), emit: bed
    script:
        """
        /usr/bin/env Rscript ${projectDir}/bin/merge_depths.R ${sampleID} ${nanosv_depths} ${sniffles_depths} ${survivor_ids} ${summary_table}
        """
}