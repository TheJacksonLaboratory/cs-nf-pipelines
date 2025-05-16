process DEEPTOOLS_PLOTHEATMAP {
        tag "$sampleID"

        cpus 8
        memory 10.GB
        time '04:00:00'

        publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
        "${params.pubdir}/${type + sampleID + '/deeptools'}"
        }, pattern: "*.pdf", mode: 'copy'

    
        container 'quay.io/biocontainers/deeptools:3.3.2--py_1'

        input:
        tuple val(sampleID), file(matrix)

        output:
        tuple val(sampleID), path("*.pdf"), emit: pdf
        tuple val(sampleID), path("*.tab"), emit: table

        script:
        """
        plotHeatmap --matrixFile $matrix \\
            --outFileName ${sampleID}.plotHeatmap.pdf \\
            --outFileNameMatrix ${sampleID}.plotHeatmap.mat.tab
        """
}
