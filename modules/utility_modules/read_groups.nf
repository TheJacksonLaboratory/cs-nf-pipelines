process READ_GROUPS {
    tag "${sampleID}"

    cpus 1
    memory 5.GB
    time '01:00:00'

    container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

    input:
        tuple val(sampleID), file(fq_reads)

    output:
        tuple val(sampleID), file("${sampleID}_read_group.txt"), emit: read_groups

    script:
        """
        /usr/bin/env python ${projectDir}/bin/read_group_from_fastq.py -o ${sampleID}_read_group.txt ${fq_reads}[0]
        """
}