process DELLY_CALL_GERMLINE {
    tag "$sampleID"
    
    cpus = 1
    memory {bam.size() < 40.GB ? 40.GB : 80.GB}
    time {bam.size() < 40.GB ? '10:00:00' : '24:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/jaxcompsci/delly-ref_data:1.1.6--refv0.2.0'
    
    input:
    tuple val(sampleID), path(bam), path(bai)
    tuple path(fasta), path(fai)
    
    output:
    tuple val(sampleID), path("${sampleID}_Delly_geno.bcf"), emit: delly_bcf

    script:
    """
    delly call \
        -q 40 \
        -x ${params.exclude_regions} \
        -s 500 \
        -o ${sampleID}_Delly.bcf \
        -g ${fasta} ${bam}

    delly call \
        -q 40 \
        -x ${params.exclude_regions} \
        -s 500 \
        -v ${sampleID}_Delly.bcf \
        -o ${sampleID}_Delly_geno.bcf \
        -g ${fasta} ${bam}

   """
}
