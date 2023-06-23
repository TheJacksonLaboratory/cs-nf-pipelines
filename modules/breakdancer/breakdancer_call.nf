process BREAKDANCER_CALL {
    tag "$sampleID"
    
    cpus = 1
    memory {bam.size() < 40.GB ? 40.GB : 80.GB}
    time {bam.size() < 40.GB ? '10:00:00' : '24:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/breakdancer:1.4.5--h25a10a7_7'
    
    input:
        tuple val(sampleID), file(bam), file(bai)
    
    output:
        tuple val(sampleID), file("${sampleID}_BreakDancer-SV"), emit: breakdancer_sv

    script:
        """
        bam2cfg.pl ${bam} > ${sampleID}_config
        breakdancer-max -r 5 -s 50 -h ${sampleID}_config > ${sampleID}_BreakDancer-SV
        """
}