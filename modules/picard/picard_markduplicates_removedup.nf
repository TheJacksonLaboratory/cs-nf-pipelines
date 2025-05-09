process PICARD_MARKDUPLICATES {
    tag "${sampleID}"

    cpus 1
    memory 50.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", mode:'copy'

    input:
        tuple val(sampleID), file(bam)

    output:
        tuple val(sampleID), file("${sampleID}.dedup.bam"), file("${sampleID}.dedup.bai"), emit: bam_and_index
        tuple val(sampleID), file("${sampleID}.dedup_metrics.txt"), emit: dedup_metrics

    script:
        String my_mem = (task.memory-1.GB).toString()
        my_mem =  my_mem[0..-4]
        """
        picard -Xmx${my_mem}G \
            MarkDuplicates \
            --MAX_RECORDS_IN_RAM 50000 \
            --INPUT ${bam} \
            --METRICS_FILE ${sampleID}.dedup_metrics.txt \
            --TMP_DIR . \
            --ASSUME_SORT_ORDER coordinate \
            --CREATE_INDEX true \
            --REMOVE_DUPLICATES true \
            --OUTPUT ${sampleID}.dedup.bam
        """
}
