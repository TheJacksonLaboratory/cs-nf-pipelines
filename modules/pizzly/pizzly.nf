process PIZZLY {

    tag "$sampleID"

    cpus 1
    memory { 10.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy 'finish'
    maxRetries 1

    container 'quay.io/biocontainers/pizzly:0.37.3--h470a237_3'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/fusions': 'pizzly' }", pattern: "*_pizzly_results.txt", mode:'copy'

    input:
        tuple val(sampleID), path(kallisto_fusions), path(kallisto_abundance)

    output:
        tuple val(sampleID), path("*_pizzly_results.txt"), emit: pizzly_fusions

    script:
    """

    max_insert_size=`python ${projectDir}/bin/rna_fusion/compute_insert_size.py ${kallisto_abundance}`

    pizzly \
    -k 31 \
    --align-score 2 \
    --insert-size "\${max_insert_size}" \
    --cache index.cache.txt \
    --gtf ${params.gtf} \
    --fasta ${params.transcript_fasta} \
    --output ${sampleID}.pizzly ${kallisto_fusions}

    pizzly_flatten_json.py ${sampleID}.pizzly.json ${sampleID}_pizzly_results.txt

    """
}


