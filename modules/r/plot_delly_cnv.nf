process PLOT_DELLY_CNV {
    tag "$sampleID"
    
    cpus 1
    memory 1.GB
    time 1.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/r-sv_cnv_annotate:4.1.1'

    publishDir "${params.pubdir}/${sampleID + '/cnv_plots'}", pattern: "*.png", mode: 'copy'

    input:
    tuple val(sampleID), path(cov), path(seg_bed)

    output:
    path('*.png'), emit: cnv_plot

    script:
    """
    Rscript ${projectDir}/bin/pta/delly_cnv_plot.r ${cov} ${seg_bed} ${sampleID}
    """
}
