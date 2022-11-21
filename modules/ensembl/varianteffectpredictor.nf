process Ensembl_Variant_Effect_Predictor {
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'
  
  container 'ensemblorg/ensembl-vep:release_93.5'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'vep' }", pattern: "*.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(normal_germline_recalibrated_vcf)
  tuple val(sampleID), file(normal_germline_recalibrated_vcf_index)

  output:
  tuple val(sampleID), file("*_vep_annotated.vcf"), emit: normal_germline_annotated_vcf

  script:
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

    """
  vep --input_file ${sampleID}_variants_raw.recalibrated.filtered.vcf \
  --output_file ${sampleID}_vep_annotated.vcf \
  --vcf \
  --species homo_sapiens \
  --assembly GRCh38 \
  --offline \
  --cache \
  --dir_cache /opt/vep/.vep \
  --cache_version 96 \
  --force \
  --no_stats \
  --everything \
  --fasta /opt/vep/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
  """
}

