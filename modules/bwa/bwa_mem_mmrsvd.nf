process BWA_MEM {
    tag "${sampleID}"

    cpus 8
    memory 250.GB
    time '72:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_6'
    stageInMode 'copy'

    input:
        tuple val(sampleID), file(fq_reads), file(read_group)
        tuple file(fasta), file(fai)
        file(bwa_index)
    output:
        tuple val(sampleID), file("${sampleID}.sam"), emit: sam
    script:
        //if (!params.fastq2)  {
        if (params.read_type == "SE")  {
            inputfq="${fq_reads[0]}"
        }
        else {
            inputfq="${fq_reads[0]} ${fq_reads[1]}"
        }
        
        """
        bwa mem -K 100000000 -R \$(cat $read_group) -t ${task.cpus} -M ${fasta} ${inputfq} > ${sampleID}.sam
        """
}
