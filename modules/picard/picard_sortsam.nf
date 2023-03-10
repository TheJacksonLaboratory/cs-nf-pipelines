process PICARD_SORTSAM {
	
    tag "${sampleID}"

    cpus 1
    memory 4.GB
    time '2:00:00'
    container 'quay.io/biocontainers/picard:2.25.1--hdfd78af_1'

    input:
        tuple val(sampleID), file(bam)

    output:
        tuple val(sampleID), file("${sampleID}.sorted.bam"), file("${sampleID}.sorted.bai"), emit: sorted_bam


    script:
        """
        picard SortSam I=${bam} O=${sampleID}.sorted.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
        """
}