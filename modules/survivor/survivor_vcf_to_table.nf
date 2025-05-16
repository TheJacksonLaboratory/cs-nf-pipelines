process SURVIVOR_VCF_TO_TABLE {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/perl:0.1.0'

    input:
        tuple val(sampleID), file(vcf)
    output:
        tuple val(sampleID), file("${sampleID}.merged.overlap.annotated.txt"), emit: annotation
    script:
        if (params.data_type == "pacbio")
            """
            /usr/bin/env bash ${projectDir}/bin/germline_sv/surv_annot.sh ${sampleID} ${vcf} pacbio
            """
        else if (params.data_type == "illumina")
            """
            /usr/bin/env bash ${projectDir}/bin/germline_sv/surv_annot.sh ${sampleID} ${vcf} illumina
            """
        else if (params.data_type == "ont")
            """
            /usr/bin/env bash ${projectDir}/bin/germline_sv/surv_annot.sh ${sampleID} ${vcf} ont
            """            
}
