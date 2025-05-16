process BISMARK_ALIGNMENT {
    tag "$sampleID"

    cpus 20
    memory 60.GB
    time 30.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bismark:0.23.1--hdfd78af_0'

    publishDir "${params.pubdir}/${sampleID + '/alignment'}", pattern: "*.bam", mode:'copy', enabled: params.keep_intermediate
    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*report.txt", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/alignment'}", pattern: "*unmapped*", mode:'copy', enabled: params.keep_intermediate
    publishDir "${params.pubdir}/${sampleID + '/alignment'}", pattern: "*ambiguous*", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(fq_reads)

    output:
    tuple val(sampleID), path("*.bam"), emit: bam
    tuple val(sampleID), path("*report.txt"), emit: report
    tuple val(sampleID), path("*ambiguous*"), emit: ambiguous_reads
    tuple val(sampleID), path("*unmapped*"), emit: unmapped_reads

    script:

    inputfq = params.read_type == 'PE' ?  "-1 ${fq_reads[0]} -2 ${fq_reads[1]}" : "-1 ${fq_reads[0]}"
    directionality = params.non_directional ? '--non_directional': ''

    aligner = params.aligner == "hisat2" ? "--hisat2" : "--bowtie2"
    seed = params.aligner == "hisat2" ? "" : "-L ${params.seedLength}"

    """
    bismark ${aligner} --parallel 4 --bam ${directionality}  ${seed} -N ${params.seedMismatch} -minins ${params.MinInsert} -maxins ${params.MaxInsert}  --unmapped --ambiguous --genome ${params.ref_fa_index} ${inputfq}
    """
}
