process PICARD_SORTSAM {
	
    tag "${sampleID}"

    cpus 1
    memory 4.GB
    time '2:00:00'
    container 'quay.io/biocontainers/picard:2.25.1--hdfd78af_1'

    input:
        tuple val(sampleID), file(bam)
        val(suffix_string)

    output:
        tuple val(sampleID), file("${sampleID}_${suffix_string}.bam"), file("${sampleID}_${suffix_string}.bai"), emit: sorted_bam


    script:
        """
        picard SortSam I=${bam} O=${sampleID}_${suffix_string}.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
        """
}