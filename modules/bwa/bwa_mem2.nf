process BWA_MEM2 {
    tag "$sampleID"

    cpus 12
    memory 65.GB
    time 48.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : '' 
        "${params.pubdir}/${type + sampleID}"
    }, pattern: "*.sam", mode: 'copy', enabled: params.keep_intermediate


    input:
    tuple val(sampleID), path(fq_reads), file(read_groups)

    output:
    tuple val(sampleID), path("*.sam"), emit: sam

    script:

    if (params.read_type == "SE"){
        inputfq="${fq_reads[0]}"
    }
    if (params.read_type == "PE"){
        inputfq="${fq_reads[0]} ${fq_reads[1]}"
    }

    score = params.bwa_min_score ? "-T ${params.bwa_min_score}" : ''
    split_hits = params.workflow == "chipseq" ? "-M" : ''
    """
    rg=\$(cat $read_groups)
    bwa-mem2 mem -R \${rg} \
    -t $task.cpus $split_hits -B ${params.mismatch_penalty} $score ${params.ref_fa_indices} $inputfq > ${sampleID}.sam
    """
}
