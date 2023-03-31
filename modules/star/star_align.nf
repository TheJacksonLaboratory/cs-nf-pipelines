process STAR_ALIGN {

    tag "$sampleID"

    cpus 12
    memory { 84.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy 'finish'
    maxRetries 1

    container 'quay.io/biocontainers/star:2.7.8a--h9ee0642_1'

    input:
        tuple val(sampleID), path(reads)
        val(args)

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
        --sjdbGTFfile ${params.gtf} \\
        ${args}
    """
}
