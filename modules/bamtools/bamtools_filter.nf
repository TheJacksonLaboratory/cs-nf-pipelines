process BAMTOOLS_FILTER {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '12:00:00'

    container 'quay.io/biocontainers/bamtools:2.5.1--h9a82719_9'

    input:
    tuple val(sampleID), file(bam)
    file(bamtools_filter_config)

    output:
    tuple val(sampleID), file("*.sorted.bam"), emit: bam

    script:
    prefix = params.read_type == 'SE' ? "${sampleID}.mLb.clN" : "${sampleID}.mLb.flT"
    """
    bamtools filter -in ${bam} -script ${bamtools_filter_config} -out ${prefix}.sorted.bam
    """
}
