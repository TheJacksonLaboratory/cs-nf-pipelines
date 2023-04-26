process PYTHON_PARSE_SURVIVOR_IDS {

    tag "$sampleID"

    cpus 1
    memory 20.GB
    time "00:30:00"

    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'

    input:
        tuple val(sampleID), path(vcf)
    output:
        tuple val(sampleID), path("${sampleID}_survivor_depths.csv"), emit: csv
    script:

    if (params.workflow == "ont")
        """
        /usr/bin/env python ${projectDir}/bin/parse_survivor_ids.py \
        -v ${vcf} \
        -o ${sampleID}_survivor_depths.csv
        """
    else error "module relies on script that currently only supports ONT"
}