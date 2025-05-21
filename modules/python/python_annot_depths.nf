process PYTHON_ANNOT_DEPTHS {
    tag "$sampleID"

    cpus 1
    memory 20.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/pysam:0.15.2--py36h02877da_7'

    input:
        tuple val(sampleID), path(vcf), path(bed)
    output:
        tuple val(sampleID), file("${sampleID}_ONT_NS_merged_variants_depths.vcf"), emit: vcf
    script:

    if (params.data_type == "ont")
        """
        /usr/bin/env python ${projectDir}/bin/germline_sv/annot_vcf_with_depths.py \
        -v ${vcf} \
        -d ${bed} \
        -o ${sampleID}_ONT_NS_merged_variants_depths.vcf
        """
    else error "module relies on script that currently only supports ONT"
}
