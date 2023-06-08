process ARIA_DOWNLOAD {

    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '10:00:00'
    errorStrategy {(task.attempt < 5) ? { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }.call() : 'finish'}
    maxRetries 6

    container 'quay.io/jaxcompsci/aria2:1.36.0'

    input:
    tuple val(sampleID), val(meta), val(read_num), val(link)

    output:
    tuple val(sampleID), val(meta), val(read_num), path("*"), emit: file

    script:

    """
    aria2c --connect-timeout=180 --retry-wait=60 --timeout=180 ${link}
    """
}
