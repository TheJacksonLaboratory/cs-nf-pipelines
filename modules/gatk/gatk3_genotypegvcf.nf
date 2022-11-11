process GATKv3_5_GENOTYPEGVCF {
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'
  
  container 'broadinstitute/gatk3:3.5-0'
  
  
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.gvcf", mode:'copy'

  input:
  tuple val(sampleID), file(gvcf)

  output:
  tuple val(sampleID), file("*.*vcf"), emit: vcf
  tuple val(sampleID), file("*.idx"), emit: idx

  script:
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  
  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar GenomeAnalysisTK.jar \
  -T GenotypeGVCFs \
  -R ${params.ref_fa} \
  --variant ${sampleID}_variants_raw.gvcf \
  -o ${sampleID}_variants_raw.vcf
  """
}
   