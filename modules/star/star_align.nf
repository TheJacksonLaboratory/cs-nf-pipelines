process STAR_ALIGN {
    tag "$sampleID"

    cpus 12
    memory 84.GB
    time 24.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/star:2.7.8a--h9ee0642_1'

    input:
        tuple val(sampleID), path(reads)
        val(args)
        path(gtf)

    output:
        tuple val(sampleID), path('*d.out.bam'), emit: bam
        tuple val(sampleID), path('*Log.final.out'), emit: log_final
        tuple val(sampleID), path('*Log.out'), emit: log_out

        tuple val(sampleID), path('*sortedByCoord.out.bam'), optional:true, emit: bam_sorted
        tuple val(sampleID), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
        tuple val(sampleID), path('*Aligned.unsort.out.bam'), optional:true, emit: bam_unsorted
        tuple val(sampleID), path('*.tab'), optional:true, emit: tab
        tuple val(sampleID), path('*.out.junction'), optional:true, emit: junction
        tuple val(sampleID), path('*.out.sam'), optional:true, emit: sam

    script:
    """
    STAR \\
        --genomeDir ${params.star_index} \\
        --readFilesIn ${reads}  \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${sampleID}_ \\
        --sjdbGTFfile ${gtf} \\
        ${args}
    """
}
