process TRIM_GALORE {
  tag "$sampleID"

  cpus 8
  memory 16.GB
  time '06:00:00'

  container 'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/trimmed_fastq' : 'trim_galore' }", pattern: "*.fq.gz", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'fastqc' }", pattern: "*_fastqc.{zip,html}", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'trim_report' }", pattern: "*trimming_report.txt", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*_fastqc.{zip,html}"), emit: trimmed_fastqc
  tuple val(sampleID), file("*.fq.gz"), emit: trimmed_fastq
  tuple val(sampleID), file("*trimming_report.txt"), emit: trim_stats

  script:
  log.info "----- Trim Galore Running on: ${sampleID} -----"

  rrbs_flag = params.workflow == "rrbs" ? '--rrbs' : ''
  paired_end = params.read_type == 'PE' ?  '--paired' : ''
  directionality = params.non_directional ? '--non_directional': ''

  """
    trim_galore --basename ${sampleID} --cores ${task.cpus} ${paired_end} ${rrbs_flag} ${directionality} --gzip --length ${params.trimLength} -q ${params.qualThreshold}  --stringency ${params.adapOverlap}  -a ${params.adaptorSeq}  --fastqc ${fq_reads}
  """
}
