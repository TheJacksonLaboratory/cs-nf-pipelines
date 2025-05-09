process DEEPTOOLS_PLOTFINGERPRINT {
    tag "${ip} vs ${control}"

    cpus 8
    memory 10.GB
    time '04:00:00'

    publishDir "${params.pubdir}/${'immuno_precip_samples/' + ip + '_vs_' + control + '/deeptools'}", pattern: "*.pdf", mode: 'copy'
    publishDir "${params.pubdir}/${'immuno_precip_samples/' + ip + '_vs_' + control + '/deeptools'}", pattern: "*.raw.txt", mode: 'copy'
    publishDir "${params.pubdir}/${'immuno_precip_samples/' + ip + '_vs_' + control + '/deeptools'}", pattern: "*.qcmetrics.txt", mode: 'copy'

    container 'quay.io/biocontainers/deeptools:3.3.2--py_1'

    input:
    tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), file(ipbam), val(control), file(controlbam), file(ipflagstat) 

    output:
    tuple val(ip), file("*.pdf"), emit : pdf          
    tuple val(ip), file("*.raw.txt"), emit : raw      
    tuple val(ip), file("*.qcmetrics.txt"), emit : qc


    script:
    extend   = (params.read_type == 'SE' && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
    """
    plotFingerprint \\
        --bamfiles ${ipbam[0]} ${controlbam[0]} \\
        --plotFile ${ip}.plotFingerprint.pdf \\
        $extend \\
        --labels $ip $control \\
        --outRawCounts ${ip}.plotFingerprint.raw.txt \\
        --outQualityMetrics ${ip}.plotFingerprint.qcmetrics.txt \\
        --skipZeros \\
        --JSDsample ${controlbam[0]} \\
        --numberOfProcessors $task.cpus \\
        --numberOfSamples $params.fingerprint_bins
    """
}
