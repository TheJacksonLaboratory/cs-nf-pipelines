process SAMTOOLS_STATS_INSERTSIZE {
    tag "$sampleID"

    cpus 1
    memory 1.GB
    time '00:10:00'

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    input:
        tuple val(sampleID), val(meta), file(bam), file(bai), val(read_ID)

    output:
        tuple val(sampleID), env read_length, env insert_size, emit: bam

    script:

    """
    samtools stats --insert-size 8000  ${bam} | grep ^SN | cut -f 2- > insert_size.txt
    read_length=`grep "maximum length" insert_size.txt | cut -d ':' -f2 | tr -d " \t\n\r"
    insert_size=`grep "insert size average" insert_size | cut -d ':' -f2 | tr -d " \t\n\r"
    """
}


