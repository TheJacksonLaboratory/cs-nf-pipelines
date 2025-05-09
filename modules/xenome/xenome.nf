process XENOME_CLASSIFY {
    tag "$sampleID"

    cpus 8
    memory 120.GB
    time 24.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/xenome:1.0.1'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.txt", mode:'copy'

    input:
    tuple val(sampleID), path(trimmed)

    output:
    tuple val(sampleID), path("${sampleID}_sorted_human*.fastq"), emit: xenome_human_fastq
    tuple val(sampleID), path("${sampleID}_sorted_mouse*.fastq"), emit: xenome_mouse_fastq
    tuple val(sampleID), path("*.txt"), emit: xenome_stats

    script:

    if (params.read_type == "SE")
        """
        xenome classify -T 8 -P ${params.xenome_prefix} --host-name mouse --graft-name human -i ${trimmed[0]} > ${sampleID}_xenome_stats.txt
        
        fastq-sort --temporary-directory ./ --id human_1.fastq > ${sampleID}_sorted_human_1.fastq
        fastq-sort --temporary-directory ./ --id mouse_1.fastq > ${sampleID}_sorted_mouse_1.fastq

        rm human_1.fastq
        rm mouse_1.fastq
        rm ambiguous_1.fastq
        rm both_1.fastq
        rm neither_1.fastq
        """

    else if (params.read_type == "PE")
        """
        xenome classify -T 8 -P ${params.xenome_prefix} --pairs --host-name mouse --graft-name human -i ${trimmed[0]} -i ${trimmed[1]} > ${sampleID}_xenome_stats.txt
        
        fastq-sort --temporary-directory ./ --id human_1.fastq > ${sampleID}_sorted_human_1.fastq
        fastq-sort --temporary-directory ./ --id mouse_1.fastq > ${sampleID}_sorted_mouse_1.fastq
        
        fastq-sort --temporary-directory ./ --id human_2.fastq > ${sampleID}_sorted_human_2.fastq
        fastq-sort --temporary-directory ./ --id mouse_2.fastq > ${sampleID}_sorted_mouse_2.fastq

        rm human_1.fastq
        rm human_2.fastq
        rm mouse_1.fastq
        rm mouse_2.fastq
        rm ambiguous_1.fastq
        rm ambiguous_2.fastq
        rm both_1.fastq
        rm both_2.fastq
        rm neither_1.fastq
        rm neither_2.fastq
        """

    else error "${params.read_type} is invalid, specify either SE or PE"
}
