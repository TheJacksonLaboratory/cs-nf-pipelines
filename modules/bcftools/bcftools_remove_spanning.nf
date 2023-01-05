process BCFTOOLS_REMOVESPANNING {
    tag "$sampleID"

    cpus = 4
    memory = 2.GB
    time = '01:00:00'

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    tuple val(sampleID), file(vcf), file(index)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:

    """
    bcftools \
    view \
    --exclude 'ALT="*"' \
    --threads ${task.cpus} \
    -o ${sampleID}_nospanning_calls.vcf \
    ${vcf}
    """

}