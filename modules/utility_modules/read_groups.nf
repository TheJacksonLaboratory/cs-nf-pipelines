process READ_GROUPS {
    tag "${sampleID}"

    cpus 1
    memory 5.GB
    time '00:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

    input:
        tuple val(sampleID), file(fq_reads)

    output:
        tuple val(sampleID), file("${sampleID}_read_group.txt"), emit: read_groups

    script:
        """
        /usr/bin/env python ${projectDir}/bin/read_group_from_fastq.py -o ${sampleID}_read_group.txt ${fq_reads}[0]
        """
}