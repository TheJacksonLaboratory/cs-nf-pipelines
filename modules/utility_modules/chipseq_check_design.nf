process CHECK_DESIGN {
    tag "$design"
    publishDir "${params.pubdir}/parsed_samplesheets", mode: 'copy'

    input:
    path(design)

    output:
    path('design_reads.csv'), emit: sample_reads
    path('design_controls.csv'), emit: study_design

    script: 
    """
    python ${projectDir}/bin/chipseq/check_design.py $design design_reads.csv design_controls.csv
    """
}

