process ABRA2 {
    tag "$sampleID"
    
    cpus = 8
    memory = 40.GB
    time = '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/abra2:2.24--h9f5acd7_1'

    input:
    tuple val(sampleID), path(bam), path(bai)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    """
    mkdir abra_temp
    java -Xmx${my_mem}G -jar /usr/local/share/abra2-2.24-1/abra2.jar --in ${bam} --out ${bam.baseName}.abra2Realign.bam --ref ${params.ref_fa} --threads ${task.cpus} --targets ${params.target_gatk} --tmpdir ./abra_temp > abra.log
    """
}
