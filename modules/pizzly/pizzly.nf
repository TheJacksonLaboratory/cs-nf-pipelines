process PIZZLY {

    tag "$sampleID"

    cpus 12
    memory { 84.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy 'finish'
    maxRetries 1

    container 'quay.io/biocontainers/pizzly:0.37.3--h470a237_3'

    input:
        tuple val(sampleID), file(kallisto_fusions)

    output:
        tuple val(sampleID), file("*pizzly.txt"), emit: pizzly_fusions

    script:
    """
    pizzly \
    -k 31 \
    --align-score 2 \
    --insert-size 400 \
    --cache index.cache.txt \
    --gtf ${params.gtf} \
    --fasta ${params.transcript_fasta} \
    --output ${sampleID}.pizzly ${kallisto_fusions}

    pizzly_flatten_json.py ${sampleID}.pizzly.json ${sampleID}.pizzly.txt

    """
}