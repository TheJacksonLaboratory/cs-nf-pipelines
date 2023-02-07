process LUMPY_SV {
  tag "$sampleID"
  
  cpus = 1
  memory = 8.GB
  time = '03:00:00'
  
  container 'quay.io/biocontainers/lumpy-sv:0.3.1--hdfd78af_3'
  
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'lumpy-sv' }", pattern:"*.vcf", mode:'copy'

  input:
  tuple val(sampleID), val(meta), file(normal_bam), file(normal_bai), val(normal_name), file(tumor_bam), file(tumor_bai), val(tumor_name)

  output:
  tuple val(sampleID), file("*_lumpy_sv.vcf"), emit: lumpy_sv_vcf

  script:
  """
  lumpyexpress \
    -B ${tumor_bam},${normal_bam} \
    -o ${sampleID}_lumpy_sv.vcf
  """
}