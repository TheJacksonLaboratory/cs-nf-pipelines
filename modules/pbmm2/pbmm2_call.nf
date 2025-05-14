process PBMM2_CALL {

    tag "$sampleID"

    cpus 8
    memory 80.GB
    time "24:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/pbmm2:1.9.0--h9ee0642_0'

    publishDir "${params.pubdir}/${sampleID + '/alignments'}", pattern: "${sampleID}.pbmm2.aligned.bam*", mode: "copy"

    input:
        tuple val(sampleID), path(fq1)
        path(pbmm2_index)
    output:
        tuple val(sampleID), file("${sampleID}.pbmm2.aligned.bam"), file("${sampleID}.pbmm2.aligned.bam.bai"), emit: pbmm2_bam 
    script:

        if (params.pbmode == "CCS")
            """
            pbmm2 align ${pbmm2_index} ${fq1} ${sampleID}.pbmm2.aligned.bam --preset CCS --sort -j ${task.cpus}
            """
        else if (params.pbmode == "CLR")
            """
            pbmm2 align ${pbmm2_index} ${fq1} ${sampleID}.pbmm2.aligned.bam --median-filter --sort -j ${task.cpus}
            """
}
