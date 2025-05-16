process PRESEQ {
    tag "$sampleID"

    cpus 4
    memory 20.GB
    time '20:00:00'
    errorStrategy 'ignore'

    container 'quay.io/biocontainers/preseq:3.1.2--h445547b_2'

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), path("*.ccurve.txt"), emit: txt
    tuple val(sampleID), path("*.log")       , emit: log

    when:
    !params.skip_preseq

    script:
    pe = params.read_type == 'SE' ? '' : '-pe'
    """
    preseq lc_extrap \\
    -output ${sampleID}.ccurve.txt \\
    -verbose \\
    -bam \\
    $pe \\
    -seed 1 \\
    $bam
    cp .command.err ${sampleID}.command.log
    """
}
