process UCSC_BEDGRAPHTOBIGWIG {
    tag "$sampleID"

    cpus 8
    memory 10.GB
    time '04:00:00'

    publishDir {
      def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
      "${params.pubdir}/${type + sampleID + '/bigwig'}"
    }, pattern: "*.bigWig", mode: 'copy'

    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1'

    input:
    tuple val(sampleID), file(bedgraph)
    file(sizes)

    output:
    tuple val(sampleID), file("*.bigWig"), emit: bigwig

    script:
    """
    bedGraphToBigWig \\
        $bedgraph \\
        $sizes \\
        ${sampleID}.bigWig


    """
}

/*
IGV steps removed, re-add if IGV is needed: 
    OUTPUT: tuple val(sampleID), file("*.igv.txt"), emit: igv_txt

    SCRIPT: find * -type f -name "*.bigWig" -exec echo -e "bigwig/"{}"\\t0,0,178" \\; > ${sampleID}.bigWig.igv.txt
*/
