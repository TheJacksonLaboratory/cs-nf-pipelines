process XENOME_CLASSIFY {
    tag "$sampleID"

    cpus 8
    memory 50.GB
    time 8.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/xenome:1.0.1'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/stats': 'xenome' }", pattern: "*.txt", mode:'copy'

    input:
    tuple val(sampleID), path(trimmed)

    output:
    tuple val(sampleID), path("human*.fastq"), emit: xenome_fastq
    tuple val(sampleID), path("mouse*.fastq"), emit: xenome_mouse_fastq
    tuple val(sampleID), path("*.txt"), emit: xenome_stats

    script:

    read_input = params.read_type == 'PE' ? "-i ${trimmed[0]} -i ${trimmed[1]}" : "-i ${trimmed[0]}"
    pairs = params.read_type == 'PE' ? "--pairs" : ""

    """
    /xenome-1.0.1-r/xenome classify -T 8 -P ${params.xenome_prefix} ${pairs} --host-name mouse --graft-name human ${read_input} > ${sampleID}_xenome_stats.txt
    """
}
