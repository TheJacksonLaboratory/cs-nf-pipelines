process BWA_MEM {
    tag "$sampleID"

    cpus 12
    memory 65.GB
    time 48.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : '' 
        "${params.pubdir}/${type + sampleID}"
    }, pattern: "*.sam", mode: 'copy', enabled: params.keep_intermediate


    input:
    tuple val(sampleID), path(fq_reads), val(index), file(read_groups)

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
    bwa mem -R \${rg} \
    -t $task.cpus $split_hits -B ${params.mismatch_penalty} $score ${params.ref_fa_indices} $inputfq > ${sampleID}_${index}.sam
    """
}

// Input: "val(index)" refers to an index value for scattered input, not the required BWA mapping index. 