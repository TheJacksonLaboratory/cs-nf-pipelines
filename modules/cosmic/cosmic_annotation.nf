process COSMIC_ANNOTATION {
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
  if (params.workflow == 'sv')
    """
    python \
    ${projectDir}/bin/sv/add_cancer_gene_census.py \
    ${params.cosmic} \
    ${vcf} \
    ${sampleID}_germline_vep_cosmic_annotated.vcf
    """
  else
    """
    ${projectDir}/bin/shared/Cosmic_Annotation_hg38.pl \
    -i1 ${params.cosmic} \
    -i2 ${vcf} > ${sampleID}_cosmic_annotation.vcf
    """
}

// cosmic for 'sv' pipeline comes from: 
// curl -H "Authorization: Basic bWlrZS5sbG95ZEBqYXgub3JnOnlSYlU4OUJ0VUotS0RjZAo=" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/cancer_gene_census.csv
// the above command provides a URL for curl download
// curl "https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v97/cancer_gene_census.csv?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1672931317&Signature=PK8YAGC%2Bh9veZqc7mIZzywkOSf0%3D" --output cancer_gene_census.csv
