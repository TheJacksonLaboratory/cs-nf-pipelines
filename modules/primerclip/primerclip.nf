process PRIMERCLIP {
    tag "$sampleID"

    cpus 1
    memory 20.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/primerclip:0.3.8--h9ee0642_1'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'primerclip' }", pattern:"*.sam", mode:'copy', enabled: params.keep_intermediate
    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/stats' : 'primerclip' }", pattern:"*primerclip_runstats.log", mode:'copy'

    input:
      tuple val(sampleID), file(sam)

    output:
      tuple val(sampleID), file("*.sam"), emit: sam
      tuple val(sampleID), file("*primerclip_runstats.log"), emit: log

    script:

    """
    primerclip ${params.masterfile} ${sam} ${sam.baseName}_primerclip.sam
    """
}
