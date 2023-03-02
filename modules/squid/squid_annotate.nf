process SQUID_ANNOTATE {

    tag "$sampleID"

    cpus 1
    memory { 10.GB * task.attempt }
    time { 5.h * task.attempt }
    errorStrategy 'finish'
    maxRetries 1

    container 'docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/fusions': 'squid' }", pattern: "*.{tsv,txt}", mode:'copy'

    input:
        tuple val(sampleID), path(txt)

    output:
        tuple val(meta), path("*annotated.txt"), emit: squid_fusions_annotated

    script:
    """
    AnnotateSQUIDOutput.py ${params.gtf} ${txt} ${sampleID}_squid_fusions_annotated.txt
    """
}