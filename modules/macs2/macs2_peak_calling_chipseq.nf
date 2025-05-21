process PEAK_CALLING_CHIPSEQ {
    tag "${ip} vs ${control}"

    cpus 2
    memory 10.GB
    time '10:00:00'

    publishDir "${params.pubdir}/${'immuno_precip_samples/' + ip + '_vs_' + control + '/macs2'}", pattern: "*_peaks.*", mode: 'copy'

    container 'quay.io/biocontainers/macs2:2.2.7.1--py39hbf8eff0_4'  

    input:
    tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), file(ipbam), val(control), file(controlbam), file(ipflagstat)
    file(peak_count_header)
    file(frip_score_header)


    output:
    tuple val(antibody), val(replicatesExist), val(multipleGroups), file("*.{narrowPeak,broadPeak}"), emit: arm_peak
    tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), val(control), file("*.{narrowPeak,broadPeak}"), emit: ip_control_peak
    tuple val(antibody), val(replicatesExist), val(multipleGroups), val(ip), val(control), emit: ip_control
    tuple val(ip), file("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(ip), file("*_peaks.gappedPeak"), emit: gapped, optional: true
    tuple val(ip), file("*_peaks.xls"), emit: xls


    script:
    broad = params.narrow_peak ? '' : "--broad --broad-cutoff ${params.broad_cutoff}"
    format = params.read_type == 'SE' ? 'BAM' : 'BAMPE'
    pileup = params.save_macs_pileup ? '-B --SPMR' : ''
    fdr = params.macs_fdr ? "--qvalue ${params.macs_fdr}" : ''
    pvalue = params.macs_pvalue ? "--pvalue ${params.macs_pvalue}" : ''
    """
    macs2 callpeak \\
    -t ${ipbam[0]} \\
    -c ${controlbam[0]} \\
    $broad \\
    -f $format \\
    -g $params.macs_gsize \\
    -n $ip \\
    $pileup \\
    $fdr \\
    $pvalue \\
    --keep-dup all
    """
}
