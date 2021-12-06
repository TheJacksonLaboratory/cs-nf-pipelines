process TRIM {
  publishDir "${params.outdir}/trimmed"

  input:
  tuple val(sampleId), file(reads)

  output:
  tuple val(sampleId), file('*R1_paired.fastq.gz'), file('*R2_paired.fastq.gz')

  script:
  """
  trimmomatic \
  PE \
  ${params.fq_path}/${reads[0]} \
  ${params.fq_path}/${reads[1]} \
  ${sampleId}_R1_paired.fastq.gz \
  ${sampleId}_R1_unpaired.fastq.gz \
  ${sampleId}_R2_paired.fastq.gz \
  ${sampleId}_R2_unpaired.fastq.gz \
  LEADING:${params.t_lead} \
  TRAILING:${params.t_trail} \
  MINLEN:${params.min_len}
  """

}
