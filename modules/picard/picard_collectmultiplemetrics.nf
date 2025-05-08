process PICARD_COLLECTMULTIPLEMETRICS {
    tag "$sampleID"

    cpus = 1
    memory = 6.GB
    time = '03:00:00'

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : '' 
        "${params.pubdir}/${type + sampleID + '/stats'}"
    }, pattern: "*.CollectMultipleMetrics.*", mode: 'copy'

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*_metrics"), emit : metrics
    tuple val(sampleID), file("*.pdf"), emit : pdf

    script:
    prefix = "${sampleID}.mLb.clN"
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    picard -Xmx${my_mem}G CollectMultipleMetrics \
    INPUT=${bam[0]} \
    OUTPUT=${prefix}.CollectMultipleMetrics \
    REFERENCE_SEQUENCE=${params.ref_fa} \
    VALIDATION_STRINGENCY=LENIENT
    """
}
