process PICARD_COLLECTRNASEQMETRICS {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.txt", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.pdf", mode:'copy'

    input:
    tuple val(sampleID), file(bam), val(strand_setting)
    val(ref_flat)
    val(ribo_intervals)


    output:
    tuple val(sampleID), file("*metrics.txt"), emit: picard_metrics

    script:

    if (strand_setting == "reverse_stranded") {
        strand_setting = "SECOND_READ_TRANSCRIPTION_STRAND"
    }

    if (strand_setting == "forward_stranded") {
        strand_setting = "FIRST_READ_TRANSCRIPTION_STRAND"
    }

    if (strand_setting == "non_stranded") {
        strand_setting = "NONE"
    }

    """
    picard CollectRnaSeqMetrics \
    I=${bam} \
    O=${sampleID}_picard_aln_metrics.txt \
    REF_FLAT=${ref_flat} \
    RIBOSOMAL_INTERVALS=${ribo_intervals} \
    STRAND=${strand_setting} \
    CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf
    """
}
