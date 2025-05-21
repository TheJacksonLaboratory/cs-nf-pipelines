process PYTHON_PARSE_DEPTHS {
    tag "$sampleID"

    cpus 1
    memory 20.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'

    input:
        tuple val(sampleID), path(vcf)
        val(caller)
    output:
        tuple val(sampleID), path("${sampleID}_${caller}_depths.csv"), emit: csv
    script:

    if (params.data_type == "ont")
        """
        /usr/bin/env python ${projectDir}/bin/germline_sv/parse_caller_depths.py \
        -v ${vcf} \
        -o ${sampleID}_${caller}_depths.csv
        """
    else error "module relies on script that currently only supports ONT"
}
