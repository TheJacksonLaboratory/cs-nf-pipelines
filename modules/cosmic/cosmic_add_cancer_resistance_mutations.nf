process COSMIC_CANCER_RESISTANCE_MUTATION {
    tag "$sampleID"

    cpus 1
    memory { 5.GB * task.attempt }
    time {1.hour * task.attempt}
    errorStrategy 'retry'
    maxRetries 1

    container 'quay.io/jaxcompsci/python-yaml:3.9.7_ps'

    input:
    tuple val(sampleID), file(vcf)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:
    """
    python \
    ${projectDir}/bin/sv/add_cancer_resistance_mutations.py \
    ${params.cosmic_cancer_resistance_muts} \
    ${vcf} \
    ${sampleID}_germline_vep_cosmic_cancerResitMut_annotated.vcf
    """
}

// cosmic for 'sv' pipeline comes from: 
// curl -H "Authorization: Basic bWlrZS5sbG95ZEBqYXgub3JnOnlSYlU4OUJ0VUotS0RjZAo=" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/CosmicResistanceMutations.tsv.gz
// the above command provides a URL for curl download
// curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v97/CosmicResistanceMutations.tsv.gz?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1672933745&Signature=nQ9AFGONT4rDKfM4UZ1cmN4J%2F%2BM%3D" --output CosmicResistanceMutations.tsv.gz