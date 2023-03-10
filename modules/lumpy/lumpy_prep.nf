process LUMPY_PREP {
    tag "$sample_name"
    
    cpus = 8
    memory = 40.GB
    time = "10:00:00"

    container 'quay.io/jaxcompsci/lumpy-ref_data:0.3.1--2'

    publishDir "${params.outdir}/alignments/mapped_lumpy", pattern: "*_alignBWA_lumpy.bam", mode: 'copy'
    
    input:
        tuple val(sampleID), file(bam), file(bai)
    
    output:
        tuple val(sampleID), file("${sampleID}_alignBWA_lumpy.bam"), emit: bam_bwa_lumpy
        tuple val(sampleID), file("${sampleID}_discordants.unsorted.bam"), emit: dis_unsorted_bam

    script:
    """
        # Clipped_rc reads mapping to Genome
        samtools sort -n ${bam} -o ${sampleID}_alignBWA_ReadNameSort.bam -@ ${params.threads}
        # manual read group info
        samtools view -h ${sampleID}_alignBWA_ReadNameSort.bam \
        | samblaster --acceptDupMarks --excludeDups --addMateTags \
                    --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20 \
        | samtools view -@ ${params.threads} -S -b - > ${sampleID}_alignBWA_lumpy.bam
        # Extract the discordant pairedend alignments
        samtools view -@ ${params.threads} -b -F 1294 ${sampleID}_alignBWA_lumpy.bam > ${sampleID}_discordants.unsorted.bam
    """
}