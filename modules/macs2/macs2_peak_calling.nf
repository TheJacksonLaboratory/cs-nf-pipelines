process PEAK_CALLING {
    tag "$sampleID"

    cpus 2
    memory 10.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/macs2:2.2.7.1--py39hbf8eff0_4'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*_peaks.narrowPeak", mode: 'copy'
    publishDir "${params.pubdir}/${sampleID}", pattern: "*_summits.bed", mode: 'copy'
    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.log", mode: 'copy'
    
    input:
    tuple val(sampleID), file(processed_bams)

    output:
    tuple val(sampleID), file("*_peaks.narrowPeak"), emit: np
    tuple val(sampleID), file("*_summits.bed")
    tuple val(sampleID), file("*_macs2.log")


    script:
    String genome = params.gen_org == 'human' ? 'hs' : 'mm'
    """
    macs2 callpeak \
    -f BAMPE \
    --nomodel \
    -g ${genome} \
    --keep-dup all \
    --cutoff-analysis \
    --tempdir ./ \
    -n ${sampleID} \
    -t ${processed_bams[0]} \
    --outdir . \
    > ${sampleID}_macs2.log 2>&1
    """
}
