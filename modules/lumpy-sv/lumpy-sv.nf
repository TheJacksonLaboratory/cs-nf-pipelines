process LUMPY_SV {
  tag "$meta.patient"
  
  cpus = 1
  memory = 8.GB
  time = '03:00:00'
  
  container 'quay.io/biocontainers/lumpy-sv:latest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'lumpy-sv' }", pattern:".vcf", mode:'copy'

  input:
  tuple val(sampleID), file(normal_bam), file(normal_bai), val(meta)
  tuple val(sampleID), file(tumor_bam), file(tumor_bai), val(meta)

  output:
  tuple val(meta), file("*.lumpy-sv.vcf"), emit: lumpy_vcf

  script:
  log.info "----- Lumpy-SV Running on: ${sampleID} -----"

  """
  lumpyexpress \
    -B ${tumor_bam},${normal_bam} \
    -o ${sampleID}.lumpy-sv.vcf
  """
}