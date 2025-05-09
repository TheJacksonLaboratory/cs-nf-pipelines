process PICARD_SORTSAM {
    tag "$sampleID"

    cpus 1
    memory { sam.size() < 60.GB ? 30.GB : 60.GB }
    time { sam.size() < 60.GB ? '10:00:00' : '24:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*_sortsam.bam", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(sam)
    val(sort_order)

    output:
    tuple val(sampleID), file("*_sortsam.bam"), emit: bam
    tuple val(sampleID), file("*_sortsam.bai"), emit: bai, optional: true

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    index_creation = sort_order == 'coordinate' ? 'CREATE_INDEX=true' : ''

    """
    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp SortSam \
    SO=${sort_order} \
    INPUT=${sam} \
    OUTPUT=${sam.baseName}_sortsam.bam \
    TMP_DIR=`pwd`/tmp \
    VALIDATION_STRINGENCY=SILENT \
    ${index_creation}
    """
}
