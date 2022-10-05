process CHECK_DESIGN {
    tag "$design"
    publishDir "${params.pubdir}/pipeline_info", mode: 'copy'

    input:
    file(design)

    output:
    file('design_reads.csv')
    file('design_controls.csv')

    script:  // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    python ${projectDir}/bin/chipseq/check_design.py $design design_reads.csv design_controls.csv
    """
}

