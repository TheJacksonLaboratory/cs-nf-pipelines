process ARRIBA {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time 2.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
  
    container 'quay.io/biocontainers/arriba:2.4.0--ha04fe3b_0'

    publishDir "${params.pubdir}/${sampleID + '/fusions' }", pattern: "*.{tsv,txt}", mode:'copy'

    input:
        tuple val(sampleID), path(bam), path(bai)
        path(gtf)

    output:
        tuple val(sampleID), path("*_arriba_fusions.tsv"), emit: arriba_fusions
        tuple val(sampleID), path("*_arriba_fusions_discarded.tsv"), emit: arriba_fusions_fail
    
    script:

    """
    arriba \\
        -x ${bam} \\
        -a ${params.fasta} \\
        -g ${gtf} \\
        -o ${sampleID}_arriba_fusions.tsv \\
        -O ${sampleID}_arriba_fusions_discarded.tsv \\
        -b ${params.arriba_blacklist} \\
        -k ${params.arriba_known_fusions} \\
        -t ${params.arriba_known_fusions} \\
        -p ${params.arriba_protein_domains}
    """
}

/*
From the documentation: 
    Note: In this execution, the same file is passed to the parameters -k and -t, because it is used for two purposes: 
    applying sensitive filtering parameters to known fusions (-k) and tagging known fusions in the tags column (-t).
    However, it is possible to use different files for these two parameters if a user wants to separate the two tasks.
*/