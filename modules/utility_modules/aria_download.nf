process ARIA_DOWNLOAD {

    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '03:00:00'

    container 'quay.io/jaxcompsci/aria2:1.36.0'

    input:
    tuple val(sampleID), val(meta), val(read_num), val(link)

    output:
    tuple val(sampleID), val(meta), val(read_num), path("*"), emit: file

    script:

    """
    aria2c ${link}
    """
}
