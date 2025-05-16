process READ_GROUPS {
    tag "$sampleID"

    cpus 1
    memory 5.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/py3_bgzip:3.10.12'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*read_group.txt", mode:'copy', enabled: params.workflow == 'rnaseq' || params.keep_intermediate

    input:
    tuple val(sampleID), path(fq_reads)
    val(picard)

    output:
    tuple val(sampleID), path("*.txt"), emit: read_groups

    script:
    if (picard=="picard"){
        p='-p'
    }
    else{
        p=''
    }
    """
    python ${projectDir}/bin/shared/read_group_from_fastq.py $p -s ${sampleID} -o ${sampleID}_read_group.txt ${fq_reads[0]}
    """
}
