process MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '04:00:00'

    container 'quay.io/biocontainers/r-base:3.5.1'

    input:
    tuple val(sampleID), file(spp), file(rdata)
    file(nsc_header)
    file(rsc_header)
    file(correlation_header)

    output:
    tuple val(sampleID), file("*.spp_nsc_mqc.tsv")        , emit: nsc
    tuple val(sampleID), file("*.spp_rsc_mqc.tsv")        , emit: rsc
    tuple val(sampleID), file("*.spp_correlation_mqc.tsv"), emit: correlation

    script:
    """
    cp $correlation_header ${sampleID}.spp_correlation_mqc.tsv
    Rscript --max-ppsize=500000 -e "load('$rdata'); write.table(crosscorr\\\$cross.correlation, file=\\"${sampleID}.spp_correlation_mqc.tsv\\", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE,append=TRUE)"

    awk -v OFS='\t' '{print "${sampleID}", \$9}'  $spp | cat $nsc_header - > ${sampleID}.spp_nsc_mqc.tsv
    awk -v OFS='\t' '{print "${sampleID}", \$10}' $spp | cat $rsc_header - > ${sampleID}.spp_rsc_mqc.tsv
    """
}
