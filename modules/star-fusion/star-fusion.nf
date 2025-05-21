process STAR_FUSION {
    tag "$sampleID"

    cpus 12
    memory 42.GB
    time 5.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'trinityctat/starfusion:1.12.0'

    publishDir "${params.pubdir}/${sampleID + '/fusions'}", pattern: "*.{tsv,txt}", mode:'copy'

    input:
        tuple val(sampleID), file(reads)

    output:
        tuple val(sampleID), file("*_star-fusion_fusions.tsv"), emit: star_fusion_fusions
        tuple val(sampleID), file("*_abridged.tsv"), emit: star_fusion_fusions_abridge
        tuple val(sampleID), file("*_abridged.coding_effect.tsv"), optional: true, emit: star_fusion_abridge_coding

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

    mv star-fusion.fusion_predictions.tsv ${sampleID}_star-fusion_fusions.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${sampleID}_star-fusion_abridged.tsv
    mv star-fusion.fusion_predictions.abridged.coding_effect.tsv ${sampleID}_star-fusion_abridged.coding_effect.tsv
    """
}

//`--readFilesCommand zcat` this option is included in STAR if files are compressed. 

/*

To build a new reference set: 

    export TMPDIR=/flashscratch/lloydm/tmp
    
    wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz --no-check-certificate
    wget https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases/download/v0.3.0/fusion_lib.Mar2021.dat.gz -O CTAT_HumanFusionLib_Mar2021.dat.gz --no-check-certificate
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/AnnotFilterRule.pm -O AnnotFilterRule.pm --no-check-certificate
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm --no-check-certificate
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3f --no-check-certificate
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3i --no-check-certificate
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3m --no-check-certificate
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3p --no-check-certificate
    gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm
    
    singularity exec /projects/omics_share/meta/containers/trinityctat-starfusion-1.12.0.img \
    /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl \
        --genome_fa /projects/compsci/omics_share/human/GRCh38/transcriptome/indices/ensembl/Homo_sapiens.GRCh38.102.all.fa \
        --gtf /projects/compsci/omics_share/human/GRCh38/transcriptome/indices/ensembl/Homo_sapiens.GRCh38.102.chr.gtf \
        --annot_filter_rule AnnotFilterRule.pm \
        --fusion_annot_lib CTAT_HumanFusionLib_Mar2021.dat.gz \
        --pfam_db Pfam-A.hmm \
        --dfam_db homo_sapiens_dfam.hmm \
        --max_readlength 150 \
        --CPU 8
*/
