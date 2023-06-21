process DELLY_CALL {
    tag "$sampleID"
    
    cpus = 1
    memory = 40.GB
    time = "10:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/jaxcompsci/delly-ref_data:0.8.3--hf3ca161_1'
    
    input:
        tuple val(sampleID), file(bam), file(bai)
        tuple file(fasta), file(fai)
    
    output:
        tuple val(sampleID), file("${sampleID}_Delly.bcf"), emit: delly_bcf

    script:
        """
        delly call \
            -q 40 \
            -x ${params.exclude_regions} \
            -s 500 \
            -o ${sampleID}_Delly.bcf \
            -g ${fasta} ${bam}
        """
}