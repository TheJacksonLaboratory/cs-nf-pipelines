process LUMPY_SV {
  tag "$sampleID"
  
  cpus = 1
  memory { normal_bam.size() < 60.GB ? 8.GB : 24.GB }
  time { normal_bam.size() < 60.GB ? '03:00:00' : '12:00:00' }
  
  container 'quay.io/biocontainers/lumpy-sv:0.3.1--hdfd78af_3'
  
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/callers' : 'lumpy-sv' }", pattern:"*.vcf", mode:'copy'

  input:
  tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name)

  output:
  tuple val(sampleID), path("*_lumpy_sv.vcf"), val(meta), val(normal_name), val(tumor_name), val('lumpy'), emit: lumpy_sv_vcf

  script:
  """
  lumpyexpress \
    -B ${tumor_bam},${normal_bam} \
    -o ${sampleID}_lumpy_sv.vcf
  """
}