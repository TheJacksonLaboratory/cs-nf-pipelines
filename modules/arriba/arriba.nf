process ARRIBA {

    tag "$sampleID"

    cpus 12
    memory { 84.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy 'finish'
    //maxRetries 1

    container 'TBD TBD TBD TBD TBD'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'star-fusion' }", pattern: "*.{tsv,txt}", mode:'copy'

    input:
        tuple val(sampleID), file(reads)

    output:
        tuple val(sampleID), file("${sampleID}_star-fusion.tsv"), optional: true, emit: star_fusion_fusions
        tuple val(sampleID), file("${sampleID}_abridged.tsv"), optional: true, emit: star_fusion_fusions_abridge
        tuple val(sampleID), file("${sampleID}_abridged.coding_effect.tsv"), optional: true, emit: star_fusion_abridge_coding

    script:

    """
    arriba \\
        -x $bam \\
        -a $fasta \\
        -g $gtf \\
        -o ${prefix}.fusions.tsv \\
        -O ${prefix}.fusions.discarded.tsv \\
        $blacklist \\
        $known_fusions \\
        $structural_variants \\
        $tags \\
        $protein_domains \\
        $args
    """
}



    withName: STAR_FOR_ARRIBA {
        publishDir = [
            path: { "${params.outdir}/star_for_arriba" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
        ]
        ext.args = '--readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outBAMcompression 0 \
        --outFilterMultimapNmax 50 \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 \
        --chimOutType WithinBAM HardClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 50'
    }
