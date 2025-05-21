process PYTHON_ANNOT_ON_TARGET {
    tag "$sampleID"

    cpus 1
    memory 20.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/pysam:0.15.2--py36h02877da_7'

    publishDir "${params.pubdir}/${sampleID}", mode:'copy'

    input:
        tuple val(sampleID), path(vcf)
    output:
        tuple val(sampleID), path("${sampleID}_ONT_NS_struct_var.vcf"), emit: vcf
     
    script:

    if (params.data_type == "ont" && params.targ_chr && params.targ_start && params.targ_end)
        """
        /usr/bin/env python ${projectDir}/bin/germline_sv/annot_vcf_with_on_target.py \
            -v ${vcf} \
            -c ${params.targ_chr} \
            -s ${params.targ_start} \
            -e ${params.targ_end} \
            -o ${sampleID}_ONT_NS_struct_var.vcf
        """
    else if (params.data_type == "ont")
        """
        /usr/bin/env python ${projectDir}/bin/germline_sv/annot_vcf_with_on_target.py \
            -v ${vcf} \
            -c NA \
            -s NA \
            -e NA \
            -o ${sampleID}_ONT_NS_struct_var.vcf
        """
    else error "module relies on script that currently only supports ONT"
}
