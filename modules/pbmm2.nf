process BUILDPBMM2INDEX {

    tag "$sampleID"

    cpus 8
    memory 250.GB
    time 72.h
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/biocontainers/pbmm2:1.9.0--h9ee0642_0'

    input:
        val(sampleID)
        file(fasta)
    output:
        file("${fasta.baseName}.mmi"), emit: pbmm2_index
    script:
        """
        pbmm2 index ${fasta} ${fasta.baseName}.mmi
        """
}

process PBMM2MAPCCS {

    tag "$sampleID"

    cpus 8
    memory 250.GB
    time 72.h
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/biocontainers/pbmm2:1.9.0--h9ee0642_0'

    input:
        val(sampleID)
        file(fq1)
        file(mmi)
    output:
        tuple val(sampleID), file(pbmm2_bam), file(pbmm2_bai), emit: pbmm2_ccs 
    script:
        """
		pbmm2 align ${mmi} ${fq1} ${$sampleID}.pbmm2.aligned.bam --preset CCS --sort -j ${task.cpus}
		"""
        
}

process PBMM2MAPCLR {
    tag "$sampleID"

    cpus 8
    memory 250.GB
    time 72.h
    maxRetries 1
    errorStrategy 'retry'

    container 'quay.io/biocontainers/pbmm2:1.9.0--h9ee0642_0'

    input:
        val(sampleID)
        file(fq1)
        file(mmi)
    output:
        tuple val(sampleID), file(pbmm2_bam), file(pbmm2_bai), emit: pbmm2_clr
    script:
        """
        pbmm2 align ${mmi} ${fq1} ${name_string}.pbmm2.aligned.bam --median-filter --sort -j ${task.cpus}
        """
}