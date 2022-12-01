process DEEPTOOLS_PLOTPROFILE {
    tag "$sampleID"

    cpus 8
    memory 10.GB
    time '04:00:00'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'deeptools' }", pattern: "*.pdf", mode: 'copy'

    container 'quay.io/biocontainers/deeptools:3.3.2--py_1'


    input:
    tuple val(sampleID), file(matrix)

    output:
    tuple val(sampleID), path("*.pdf"), emit: pdf
    tuple val(sampleID), path("*.tab"), emit: table

    script:
    """
    plotProfile --matrixFile ${sampleID}.computeMatrix.mat.gz \\
        --outFileName ${sampleID}.plotProfile.pdf \\
        --outFileNameData ${sampleID}.plotProfile.tab
    """
}
