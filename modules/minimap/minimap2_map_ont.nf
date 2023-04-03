process MINIMAP2_MAP_ONT {
    tag "$sampleID"
    cpus 16
    memory 24.GB
    time "24:00:00"

    container 'quay.io-biocontainers-minimap2-2.24--h7132678_1'

    publishDir "${params.pubdir}/", pattern: "*_aln.sam", mode:'copy', enabled: params.keep_intermediate

    input:
        tuple val(sampleID), file(fastq)
        path(minimap2_index)
    output:
        tuple val(sampleID), file("${sampleID}_aln.sam"), emit: minimap_sam
    script:
        """
        minimap2 --secondary=no -t ${task.cpus} --MD -ax map-ont ${minimap2_index} ${fastq}  > ${sampleID}_aln.sam
        """
29
}