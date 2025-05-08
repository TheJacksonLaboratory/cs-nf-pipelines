process PIZZLY {

    tag "$sampleID"

    cpus 1
    memory 10.GB
    time 2.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/pizzly:0.37.3--h470a237_3'

    publishDir "${params.pubdir}/${sampleID + '/fusions'}", pattern: "*_pizzly_fusions.txt", mode:'copy'

    input:
        tuple val(sampleID), path(kallisto_fusions), path(kallisto_insert_size)
        path(gtf)

    output:
        tuple val(sampleID), path("*_pizzly_fusions.txt"), emit: pizzly_fusions

    script:
    """
    insert_size="\$(cat ${kallisto_insert_size})"

    pizzly \
    -k 31 \
    --align-score 2 \
    --insert-size "\${insert_size}" \
    --cache index.cache.txt \
    --gtf ${gtf} \
    --fasta ${params.transcript_fasta} \
    --output ${sampleID}.pizzly ${kallisto_fusions}

    pizzly_flatten_json.py ${sampleID}.pizzly.json ${sampleID}_pizzly_fusions.txt
    """
}
