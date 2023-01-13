process GATK_DEPTHOFCOVERAGE {

  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '05:00:00'

  container 'broadinstitute/gatk:4.2.4.1'
  file(params.ref_fai)

  input:
  tuple val(sampleID), file(bam), file(bai)
  val(L)

  output:
  tuple val(sampleID), file("*_gatk_temp.txt"), emit: txt

  script:
  log.info "----- GATK Depth of Coverage Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" DepthOfCoverage \
  -R ${params.ref_fa} \
  --output-format TABLE \
  -O ${sampleID}_gatk_temp.txt \
  -I ${bam} \
  -L  ${L} \
  --omit-per-sample-statistics \
  --omit-interval-statistics \
  --omit-locus-table \
  """
}