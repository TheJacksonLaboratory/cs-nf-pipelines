process GATKv3_5_ApplyRecalibration {
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'
  
  container 'broadinstitute/gatk3:3.5-0'
  
  
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.*recalibrated.filtered.vcf"), emit: vcf


  script:
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  
  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -R ${params.ref_fa} \
  -input ${sampleID}_variants_raw.vcf \
  --ts_filter_level 99.6 \
  -tranchesFile ${sampleID}.tranches \
  -recalFile ${sampleID}.recal \
  -mode SNP
  -o ${sampleID}_variants_raw.recalibrated.filtered.vcf
  """
}