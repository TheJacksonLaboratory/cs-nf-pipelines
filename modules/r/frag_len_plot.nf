process FRAG_LEN_PLOT {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*fraglen_plot.pdf", mode: 'copy'
    container 'quay.io/jaxcompsci/rstudio:4.2.0' 

    input:
    tuple val(sampleID), file(frag_len_count)

    output:
    tuple val(sampleID), file("*fraglen_plot.pdf")
    tuple val(sampleID), file("*_spline_table.txt"), emit: spline_table

    script:
    """
    Rscript ${projectDir}/bin/atac/fragment_length_plot.R ${frag_len_count} ${sampleID}_spline_table.txt
    mv fraglen_plot.pdf ${sampleID}_fraglen_plot.pdf
    """
}
