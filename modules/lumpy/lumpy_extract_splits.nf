process LUMPY_EXTRACT_SPLITS {
    tag "$sample_name"
    
    cpus = 8
    memory = 40.GB
    time = "10:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/jaxcompsci/lumpy-ref_data:0.3.1--refv0.2.0'
    
    input:
        tuple val(sampleID), path(bam)
    
    output:
        tuple val(sampleID), path("${sampleID}_splitreads.bam"), emit: bam_bwa_lumpy

    script:
    """
        samtools view -h ${bam} \
        | extractSplitReads_BwaMem -i stdin \
        | samtools view -Sb - > ${sampleID}_splitreads.bam
    """
}
