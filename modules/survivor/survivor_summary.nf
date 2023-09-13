process SURVIVOR_SUMMARY {

    tag "$sampleID"

    cpus 1
    memory 2.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/biopython-pyvcf:1.78'

    input:
        tuple val(sampleID), file(vcf)
    output:
        tuple val(sampleID), file("${sampleID}.survivor_summary.csv"), emit: csv
    script:
        """
        /usr/bin/env python ${projectDir}/bin/germline_sv/sv_to_table.py -v ${vcf} -o ${sampleID}.survivor_summary.csv
        """
}
