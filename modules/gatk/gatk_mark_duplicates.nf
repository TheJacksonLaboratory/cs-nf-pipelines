process GATK_MARK_DUPLICATES {
    tag "${sampleID}"

    cpus 1
    memory 50.GB
    time '10:00:00'
    container 'quay.io/biocontainers/gatk4:4.1.8.1--py38_0'

    publishDir "${params.pubdir}/alignments", mode:'copy'

    input:
        tuple val(sampleID), file(bam)

    output:
        tuple val(sampleID), file("${sampleID}.md.bam"), file("${sampleID}.md.bai"), emit: bam_and_index
        tuple val(sampleID), file("${sampleID}.md.metrics"), emit: dedup_metrics

    script:
        String my_mem = (task.memory-1.GB).toString()
        my_mem =  my_mem[0..-4]
        """
        gatk --java-options -Xmx${my_mem}G \
            MarkDuplicates \
            --MAX_RECORDS_IN_RAM 50000 \
            --INPUT ${bam} \
            --METRICS_FILE ${sampleID}.md.metrics \
            --TMP_DIR . \
            --ASSUME_SORT_ORDER coordinate \
            --CREATE_INDEX true \
            --REMOVE_DUPLICATES true \
            --OUTPUT ${sampleID}.md.bam
        """
}