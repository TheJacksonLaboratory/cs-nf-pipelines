process PICARD_CLEANSAM {
    tag "$sampleID"

    cpus = 1
    memory { bam.size() < 60.GB ? 10.GB : 30.GB }
    time { bam.size() < 60.GB ? '10:00:00' : '20:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*.cleaned.bam", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*.cleaned.bam"), emit: cleaned_bam

    script:

    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    picard -Xmx${my_mem}G CleanSam \
    I=${bam} \
    TMP_DIR=${workDir}/temp \
    VALIDATION_STRINGENCY=SILENT \
    O=${bam.baseName}.cleaned.bam
    """
}
