process PBSV_DISCOVER {

    tag "$sampleID"

    cpus 8
    memory 40.GB
    time "4:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/pbsv-td_refs:2.8.0--refv0.2.0'

    input:
        tuple val(sampleID), file(pbmm2_bam), file(pbmm2_bai)
    output:
        tuple val(sampleID), file("${sampleID}.svsig.gz"), emit: pbsv_svsig
    script:
        if (params.pbsv_tandem)
            """
            pbsv discover --tandem-repeats ${params.tandem_repeats} ${pbmm2_bam} ${sampleID}.svsig.gz
            """
        else
            """
             pbsv discover ${pbmm2_bam} ${sampleID}.svsig.gz
            """
}
