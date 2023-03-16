process GATK_HAPLOTYPECALLER_INTERVAL {
// Note about this module
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time = '11:00:00'

  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(sampleID), file(bam), file(bai), val(chrom)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf
  tuple val(sampleID), file("*.idx"), emit: idx

  script:

  log.info "----- GATK Haplotype Caller Running on Chromosome ${chrom} for sample: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" HaplotypeCaller  \
  -R ${params.ref_fa} \
  -I ${bam} \
  -O ${sampleID}_HaplotypeCaller_${chrom}.vcf \
  -L ${chrom} \
  -stand-call-conf ${params.call_val}
  """
}