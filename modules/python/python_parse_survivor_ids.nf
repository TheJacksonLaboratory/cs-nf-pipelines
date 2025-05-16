process PYTHON_PARSE_SURVIVOR_IDS {
    tag "$sampleID"

    cpus 1
    memory 20.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'

    input:
        tuple val(sampleID), path(vcf)
    output:
        tuple val(sampleID), path("${sampleID}_survivor_depths.csv"), emit: csv
    script:

    if (params.data_type == "ont")
        """
        /usr/bin/env python ${projectDir}/bin/germline_sv/parse_survivor_ids.py \
        -v ${vcf} \
        -o ${sampleID}_survivor_depths.csv
        """
    else error "module relies on script that currently only supports ONT"
}
