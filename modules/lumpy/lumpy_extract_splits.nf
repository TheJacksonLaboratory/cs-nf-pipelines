process LUMPY_EXTRACT_SPLITS {
    tag "$sample_name"
    
    cpus = 8
    memory = 40.GB
    time = "10:00:00"

    container 'quay.io/jaxcompsci/lumpy-ref_data:0.3.1--2'

    publishDir "${params.outdir}/alignments/mapped_lumpy", pattern: "*_alignBWA_lumpy.bam", mode: 'copy'
    
    input:
        tuple val(sampleID), file(bam), file(bai)
    
    output:
        tuple val(sampleID), file("${sampleID}_splitreads.bam"), emit: bam_bwa_lumpy

    script:
    """
        samtools view -h ${bam} \
        | extractSplitReads_BwaMem -i stdin \
        | samtools view -Sb - > ${sampleID}_splitreads.bam
    """
}