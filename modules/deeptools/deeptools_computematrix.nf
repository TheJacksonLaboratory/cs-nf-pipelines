process DEEPTOOLS_COMPUTEMATRIX {
    tag "$sampleID"

    cpus 8
    memory 10.GB
    time '04:00:00'

    container 'quay.io/biocontainers/deeptools:3.3.2--py_1'

    input:
    tuple val(sampleID), file(bigwig)
    file(bed)

    output:
    tuple val(sampleID), file("*.mat.gz") , emit: matrix
    tuple val(sampleID), file("*.mat.tab"), emit: table

    script:
    """
    computeMatrix scale-regions \\
        --regionsFileName $bed \\
        --scoreFileName $bigwig \\
        --outFileName ${sampleID}.computeMatrix.mat.gz \\
        --outFileNameMatrix ${sampleID}.computeMatrix.vals.mat.tab \\
        --regionBodyLength 1000 \\
        --beforeRegionStartLength 3000 \\
        --afterRegionStartLength 3000 \\
        --skipZeros \\
        --smartLabels \\
        --numberOfProcessors $task.cpus
    """
}
