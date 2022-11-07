process GATK_SELECTVARIANTS {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '03:00:00'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(vcf), file(idx)
  val(indel_snp)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf
  tuple val(sampleID), file("*.idx"), emit: idx

  script:
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" SelectVariants \
  -R ${params.ref_fa} \
  -V ${vcf} \
  -select-type ${indel_snp} \
  -O ${sampleID}_selectedvariants_${indel_snp}.vcf
  """
}