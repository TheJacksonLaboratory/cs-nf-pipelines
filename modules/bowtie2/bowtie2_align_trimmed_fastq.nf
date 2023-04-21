process ALIGN_TRIMMED_FASTQ {
  tag "$sampleID"

  cpus 16
  memory 30.GB
  time '48:00:00'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'bowtie2' }", pattern: "*.log", mode: 'copy'
  container 'biocontainers/bowtie2:v2.4.1_cv1'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*.sam"), emit: sam
  tuple val(sampleID), file("*_bowtie2.log"), emit: bowtie_log

  script:
  String options = params.bowtieVSensitive  == 'true' ? '--very-sensitive' : ''
  """
  bowtie2 \
  ${options} \
  -X ${params.bowtieMaxInsert} \
  -q \
  -p $task.cpus \
  -x ${params.bowtie2Index} \
  -1 ${fq_reads[0]} \
  -2 ${fq_reads[1]} \
  -S ${sampleID}.sam \
  2>${sampleID}_bowtie2.log \
  """

}
