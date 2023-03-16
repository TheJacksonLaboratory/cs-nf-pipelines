process SURVIVOR_VCF_TO_TABLE {

    tag "$sampleID"

    cpus 1
    memory 2.GB
    time "00:30:00"

    container 'quay.io/jaxcompsci/perl:0.1.0'

    input:
        tuple val(sampleID), file(vcf)
    output:
        tuple val(sampleID), file("${sampleID}.merged.overlap.annotated.txt"), emit: annotation
    script:
        if (params.workflow == "pacbio_ccs")
            """
            /usr/bin/env bash ${projectDir}/bin/surv_annot.sh ${sampleID} ${vcf} pacbio
            """
        else if (params.workflow == "illumina")
            """
            /usr/bin/env bash ${projectDir}/bin/surv_annot.sh ${sampleID} ${vcf} illumina
            """
}