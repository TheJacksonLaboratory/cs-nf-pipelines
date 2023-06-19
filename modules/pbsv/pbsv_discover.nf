process PBSV_DISCOVER {

    tag "$sampleID"

    cpus 8
    memory 40.GB
    time "4:00:00"

    container 'quay.io/jaxcompsci/pbsv-td_refs:2.8.0--h9ee0642_0'

    input:
        tuple val(sampleID), file(pbmm2_bam), file(pbmm2_bai)
    output:
        tuple val(sampleID), file("${sampleID}.svsig.gz"), emit: pbsv_svsig
    script:
        if (params.pbsv_tandem)
            """
            pbsv discover --tandem-repeats ${params.tandem_repeats} ${pbmm2_bam} ${sampleID}.svsig.gz
            """
        else
            """
             pbsv discover ${pbmm2_bam} ${sampleID}.svsig.gz
            """
}