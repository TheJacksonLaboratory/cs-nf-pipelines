process CHECK_DESIGN {
    tag "$design"

    cpus 1
    memory 15.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'
    
    publishDir "${params.pubdir}/parsed_samplesheets", mode: 'copy'

    input:
    path(design)

    output:
    path('design_reads.csv'), emit: sample_reads
    path('design_controls.csv'), emit: study_design

    script: 
    """
    python ${projectDir}/bin/chipseq/check_design.py $design design_reads.csv design_controls.csv
    """
}
