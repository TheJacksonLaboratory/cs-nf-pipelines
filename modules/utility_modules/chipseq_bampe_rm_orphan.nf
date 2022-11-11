process BAMPE_RM_ORPHAN {
    tag "$sampleID"

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

