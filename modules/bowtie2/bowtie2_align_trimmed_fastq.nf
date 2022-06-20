process ALIGN_TRIMMED_FASTQ {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bowtie2' }", pattern: "*.log", mode: 'copy'
  container 'docker://biocontainers/bowtie2:v2.4.1_cv1'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*.sam"), emit: sam
  tuple val(sampleID), file("*_bowtie2.log"), emit: bowtie_log

  script:
  log.info "----- Bowtie2 Running on: ${sampleID} -----"
  """
  bowtie2 \
  --very-sensitive \
  -X ${params.bowtieMaxInsert} \
  -q \
  -p $task.cpus \
  -x ${params.bowtieIndex} \
  -1 ${fq_reads[0]} \
  -2 ${fq_reads[1]} \
  -S ${sampleID}.sam \
  2>${sampleID}_bowtie2.log \
  """

}
