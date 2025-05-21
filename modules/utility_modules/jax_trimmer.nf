process JAX_TRIMMER {
    tag "$sampleID"

    cpus 1
    memory 30.GB
    time '30:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*_stat", mode:'copy'

    input:
    tuple val(sampleID), path(fq_reads)

    output:
    tuple val(sampleID), path("*_stat"), emit: quality_stats
    tuple val(sampleID), path("*.trimmed.R*.fastq"), emit: trimmed_fastq

    script:

    if (params.read_type == "SE"){
        mode_HQ="-S -M"
        inputfq="${fq_reads[0]}"
        rename="mv ${fq_reads[0]}_filtered_trimmed ${sampleID}.trimmed.R1.fastq"
    }
    if (params.read_type == "PE"){
        mode_HQ="-M"
        inputfq="${fq_reads[0]} ${fq_reads[1]}"
        rename="mv ${fq_reads[0]}_filtered_trimmed ${sampleID}.trimmed.R1.fastq && mv ${fq_reads[1]}_filtered_trimmed ${sampleID}.trimmed.R2.fastq"
    }

    """
    python ${projectDir}/bin/shared/filter_trim.py $mode_HQ ${params.min_pct_hq_reads} -p ${params.hq_pct} -f ${params.filter_hq} -t ${params.trim_hq} $inputfq 
    ${rename}
    """
}
