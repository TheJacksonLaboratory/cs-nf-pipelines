process SNV_TO_MNV_FINAL_FILTER {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '04:00:00'

    container 'quay.io/jaxcompsci/bedtools-python3:2.26.0'

    input:
    tuple val(sampleID), file(vcf), val(meta), val(chrom)

    output:
    tuple val(sampleID), file("*.vcf"), val(meta), val(chrom), emit: vcf

    script:
    """
    python \
    ${projectDir}/bin/sv/SNVsToMNVs_CountsBasedFilter_AnnotateHighConf.py \
    -i ${vcf} \
    -o ${sampleID}_mnv_final_filtered_${chrom}.vcf
    """
}
