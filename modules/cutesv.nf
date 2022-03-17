// ccs and clr params handled in config

process cutesv_css {
    tag "$sampleID"

    cpus = 8
    time= '72:00:00'
    memory = '250.GB'
    maxRetries = 1
    clusterOptions = '-q batch'
    errorStrategy = 'retry'

    container 'quay.io/biocontainers/cutesv:1.0.9--py_0'

    input:
    tuple val(sampleID), file(bam)
    tuple val(sampleID), file(bai)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:
    log.info "----- Cute SV Running on: ${sampleID} -----"

    """
    cuteSV ${bam} ${params.fasta} ${sampleID}_cutesv_calls.vcf ./ \
    --threads ${task.cpus} \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 \
    --max_cluster_bias_DEL 1000 \
    --diff_ratio_merging_DEL 0.5
    """
}
