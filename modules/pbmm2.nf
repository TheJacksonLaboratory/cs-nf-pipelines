process BUILDPBMM2INDEX {

    tag "$sampleID"

    cpus 8
    memory { 40.GB * task.attempt }
    time { 2.h * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/biocontainers/pbmm2:1.9.0--h9ee0642_0'

    input:
        val(sampleID)
        path(fasta)
    output:
        path "${fasta.baseName}.mmi", emit: pbmm2_index
    script:
        """
        pbmm2 index ${fasta} ${fasta.baseName}.mmi
        """
}

process PBMM2MAPCCS {

    tag "$sampleID"

    cpus 8
    memory { 40.GB * task.attempt }
    time { 8.h * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/biocontainers/pbmm2:1.9.0--h9ee0642_0'

    input:
        val(sampleID)
        path(fq1)
        path(pbmm2_index)
    output:
        tuple val(sampleID), file("${sampleID}.pbmm2.aligned.bam"), file("${sampleID}.pbmm2.aligned.bam.bai"), emit: pbmm2_ccs 
    script:
        """
		pbmm2 align ${pbmm2_index} ${fq1} ${sampleID}.pbmm2.aligned.bam --preset CCS --sort -j ${task.cpus}
		"""
        
}

process PBMM2MAPCLR {
    tag "$sampleID"

    cpus 8
    memory { 40.GB * task.attempt }
    time { 8.h * task.attempt }
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/biocontainers/pbmm2:1.9.0--h9ee0642_0'

    input:
        val(sampleID)
        path(fq1)
        path(pbmm2_index)
    output:
        tuple val(sampleID), file("${sampleID}.pbmm2.aligned.bam"), file("${sampleID}.pbmm2.aligned.bam.bai"), emit: pbmm2_clr
    script:
        """
        pbmm2 align ${mmi} ${fq1} ${sampleID}.pbmm2.aligned.bam --median-filter --sort -j ${task.cpus}
        """
}