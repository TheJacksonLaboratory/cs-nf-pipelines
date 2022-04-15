// ccs and clr params handled in config

process CURESV {
    // SET THE PARAMS IN CONFIG
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
    --max_cluster_bias_INS ${params.max_cluster_bias_INS} \
    --diff_ratio_merging_INS ${params.diff_ratio_merging_INS} \
    --max_cluster_bias_DEL ${params.max_cluster_bias_DEL} \
    --diff_ratio_merging_DEL ${params.diff_ratio_merging_DEL}
    """
}
