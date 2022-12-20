process SAMTOOLS_STATS_INSERTSIZE {
    tag "$sampleID"

    cpus 1
    memory 1.GB
    time '00:10:00'

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'samtools' }", pattern: "*insert_size.txt", mode:'copy'

    input:
        tuple val(sampleID), val(meta), file(bam), file(bai), val(read_ID)

    output:
        tuple val(sampleID), env(read_length), env(insert_size), emit: read_length_insert_size
        file("*insert_size.txt")

    script:
    """
    samtools stats --insert-size 8000  ${bam} | grep ^SN | cut -f 2- > ${sampleID}_insert_size.txt
    read_length=`grep "maximum length" ${sampleID}_insert_size.txt | cut -d ':' -f2 | tr -d " \\t\\n\\r"`
    insert_size=`grep "insert size average" ${sampleID}_insert_size.txt | cut -d ':' -f2 | tr -d " \\t\\n\\r"`
    """
}


