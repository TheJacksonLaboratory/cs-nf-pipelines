process PLOT_MACS2_QC {
    
    cpus 2
    memory 10.GB
    time '10:00:00'

    container 'quay.io/biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0'
    
    publishDir "${params.pubdir}/${'immuno_precip_samples/cross_sample_plots'}", mode: 'copy'

    input:
    file(peaks)

    when:
    params.macs_gsize && !params.skip_peak_annotation && !params.skip_peak_qc

    output:
    path '*.txt'       , emit: txt
    path '*.pdf'       , emit: pdf

    script: // This script was bundled withing the nf-core/chipseq/bin/ directory
    def peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    """
    ${projectDir}/bin/chipseq/plot_macs_qc.r \\
        -i ${peaks.join(',')} \\
        -s ${peaks.join(',').replaceAll("_peaks.${peak_type}","")} \\
        -o ./ \\
        -p macs_peak
    """
}
