process COMPRESS_INDEX_MERGED_VCF {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

  input:
  tuple val(sampleID), file(vcf), val(meta)

  output:
  tuple val(sampleID), file("*.vcf.gz"), file("*.vcf.gz.tbi"), val(meta), val(normal_name), val(tumor_name), emit: compressed_vcf_tbi

  script:
    normal = meta.normal_id
    tumor = meta.tumor_id

    """
    bgzip \
    -c \
    ${vcf} \
    > ${vcf}.gz

    tabix ${vcf}.gz
    """
}
