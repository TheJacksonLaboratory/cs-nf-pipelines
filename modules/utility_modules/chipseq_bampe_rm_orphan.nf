process BAMPE_RM_ORPHAN {
    tag "$sampleID"

    container 'quay.io/biocontainers/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0'


    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam

    script:  // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    prefix = "${sampleID}.mLb.clN"
    """
    python ${projectDir}/bin/chipseq/bampe_rm_orphan.py ${bam[0]} ${prefix}.bam --only_fr_pairs
    """
}

