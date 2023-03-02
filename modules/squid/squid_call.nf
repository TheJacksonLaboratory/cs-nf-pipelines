process SQUID {

    tag "$sampleID"

    cpus 1
    memory { 10.GB * task.attempt }
    time { 5.h * task.attempt }
    errorStrategy 'finish'
    maxRetries 1

    container 'docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a'

    input:
        tuple val(sampleID), path(bam), path(chimeric_bam)

    output:
        tuple val(sampleID), path("*sv.txt"), emit: squid_fusions

    script:
    """
    squid -b ${bam} -c ${chimeric_bam} -o ${sampleID}.squid.fusions
    """
}