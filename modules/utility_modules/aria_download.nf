process ARIA_DOWNLOAD {

    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '03:00:00'

    container 'quay.io/jaxcompsci/aria2:1.36.0'

    input:
    tuple val(sampleID), val(meta), val(links)

    output:
    tuple val(sampleID), val(meta), file("*FAS*"), emit: fastq, optional: true
    // tuple val(sampleID), val(meta), file("*fas*"), emit: fastq, optional: true

    script:
    link_list = links.collect { "$it" }.join('\n')

    """
    echo "${link_list}" > link_list.txt
    aria2c -i link_list.txt
    """
}
