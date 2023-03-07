process SNPSIFT_ANNOTATE {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  // SNPEFF and SNPSIFT need updating
  container 'quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1'

  input:
  tuple val(sampleID), file(vcf)
  path(annot_source)
  path(annot_index)
  val(output_suffix)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:

  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  java -Xmx${my_mem}G -jar /opt/snpEff/SnpSift.jar \
  annotate -noDownload -id ${annot_source} ${vcf} > ${vcf.baseName}_${output_suffix}.vcf
  """
}