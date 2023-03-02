process ARRIBA {

    tag "$sampleID"

    cpus 1
    memory { 10.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy 'finish'
    // maxRetries 1

    container 'quay.io/biocontainers/arriba:2.4.0--ha04fe3b_0'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/fusions' : 'arriba' }", pattern: "*.{tsv,txt}", mode:'copy'

    input:
        tuple val(sampleID), path(bam), path(bai)

    output:
        tuple val(sampleID), path("*_arriba_results.tsv"), emit: arriba_fusions
        tuple val(sampleID), path("*_arriba_fusions_discarded.tsv"), emit: arriba_fusions_fail
    
    script:

    """
    arriba \\
        -x ${bam} \\
        -a ${params.fasta} \\
        -g ${params.gtf} \\
        -o ${sampleID}_arriba_results.tsv \\
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