process STAR_FUSION {

    tag "$sampleID"

    cpus 12
    memory { 84.GB * task.attempt }
    time { 24.h * task.attempt }
    errorStrategy 'finish'
    //maxRetries 1

    container 'trinityctat/starfusion:1.12.0'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'star-fusion' }", pattern: "*.{tsv,txt}", mode:'copy'

    input:
        tuple val(sampleID), file(reads)

    output:
        tuple val(sampleID), file("${sampleID}_star-fusion.tsv"), optional: true, emit: star_fusion_fusions
        tuple val(sampleID), file("${sampleID}_abridged.tsv"), optional: true, emit: star_fusion_fusions_abridge
        tuple val(sampleID), file("${sampleID}_abridged.coding_effect.tsv"), optional: true, emit: star_fusion_abridge_coding

    script:
    def avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    option = params.read_type == 'PE' ? "--left_fq ${reads[0]} --right_fq ${reads[1]}" : "--left_fq ${reads[0]}"
    def extra_params = params.star_fusion_opt ? params.star_fusion_opt : ''

    """
    STAR \\
        --genomeDir ${params.star_index} \\
        --readFilesIn ${reads} \\
        --twopassMode Basic \\
        --outReadsUnmapped None \\
        --chimSegmentMin 12 \\
        --chimJunctionOverhangMin 12 \\
        --alignSJDBoverhangMin 10 \\
        --alignMatesGapMax 100000 \\
        --alignIntronMax 6000 \\
        --chimSegmentReadGapMax 3 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --runThreadN ${task.cpus} \\
        --outSAMstrandField intronMotif ${avail_mem} \\
        --outSAMunmapped Within \\
        --outSAMtype BAM Unsorted \\
        --outSAMattrRGline ID:GRPundef \\
        --chimMultimapScoreRange 10 \\
        --chimMultimapNmax 10 \\
        --chimNonchimScoreDropMin 10 \\
        --peOverlapNbasesMin 12 \\
        --peOverlapMMp 0.1 \\
        --sjdbOverhang ${params.read_length - 1} \\
        --chimOutJunctionFormat 1

    STAR-Fusion \\
        --genome_lib_dir ${params.star_fusion_ref} \\
        -J Chimeric.out.junction \\
        ${option} \\
        --CPU ${task.cpus} \\
        --examine_coding_effect \\
        --output_dir . ${extra_params}

    mv star-fusion.fusion_predictions.tsv ${sampleID}_star-fusion.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${sampleID}_abridged.tsv
    mv star-fusion.fusion_predictions.abridged.coding_effect.tsv ${sampleID}_abridged.coding_effect.tsv
    """
}

//`--readFilesCommand zcat` this option is included in STAR if files are compressed. 