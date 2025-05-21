process FEATURE_COUNTS {
    tag "$sampleID"

    cpus 4
    memory 4.GB
    time '10:00:00'
    errorStrategy 'ignore'
    
    container 'quay.io/biocontainers/subread:1.6.4--h84994c4_1'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*_peaks_countMatrix.txt", mode: 'copy'
    
    input:
    tuple val(sampleID), file(processed_bams), file(peak_cvg_saf)

    output:
    tuple val(sampleID), file("*_peaks_countMatrix.txt")

    script:
    """
    featureCounts \
    -a ${peak_cvg_saf} \
    -F SAF -p \
    -T $task.cpus \
    -o ${sampleID}_peaks_countMatrix.txt \
    ${processed_bams[0]}
    """
}
