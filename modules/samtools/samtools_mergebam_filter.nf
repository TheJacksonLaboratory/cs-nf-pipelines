process SAMTOOLS_MERGEBAM_FILTER {
    tag "$sampleID"

    cpus 2
    memory 4.GB
    time '10:00:00'

    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    input:
    tuple val(sampleID), file(in_file)
    file(bed)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam

    script:
    // Setup for chipseq pipeline

    prefix = params.read_type == 'SE' ? "${sampleID}.mLb.clN" : "${sampleID}.mLb.flT"
    filter_params = params.read_type == 'SE' ? '-F 0x004' : '-F 0x004 -F 0x0008 -f 0x001'
    dup_params = params.keep_dups ? '' : '-F 0x0400'
    multimap_params = params.keep_multi_map ? '' : '-q 1'
    blacklist_params = params.blacklist ? "-L $bed" : ''

    """
    samtools view \\
        $filter_params \\
        $dup_params \\
        $multimap_params \\
        $blacklist_params \\
        -b ${in_file} > ${prefix}.bam
    """
}
