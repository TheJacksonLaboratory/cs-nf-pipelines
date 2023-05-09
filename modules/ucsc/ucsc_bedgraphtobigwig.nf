process UCSC_BEDGRAPHTOBIGWIG {
    tag "$sampleID"

    cpus 8
    memory 10.GB
    time '04:00:00'

    publishDir {
      def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
      "${params.pubdir}/${ params.organize_by=='sample' ? type+sampleID+'/bigwig' : 'ucsc'}"
    }, pattern: "*.bigWig", mode: 'copy'

    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1'

    input:
    tuple val(sampleID), file(bedgraph)
    file(sizes)

    output:
    tuple val(sampleID), file("*.bigWig"), emit: bigwig
    tuple val(sampleID), file("*.igv.txt"), emit: igv_txt

    script:
    """
    bedGraphToBigWig \\
        $bedgraph \\
        $sizes \\
        ${sampleID}.bigWig

    find * -type f -name "*.bigWig" -exec echo -e "bigwig/"{}"\\t0,0,178" \\; > ${sampleID}.bigWig.igv.txt

    """
}
