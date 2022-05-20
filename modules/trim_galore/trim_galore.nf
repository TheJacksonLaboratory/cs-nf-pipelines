process TRIM_GALORE {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'

  container 'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/trimmed_fastq' : 'trim_galore' }", pattern: "*.fq.gz", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'fastqc' }", pattern: "*_fastqc.{zip,html}", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'trim_report' }", pattern: "*trimming_report.txt}", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*_fastqc.{zip,html}"), emit: trimmed_fastqc // WHAT STATS DOES IT EMIT? 
  tuple val(sampleID), file("*.fq.gz"), emit: trimmed_fastq
  tuple val(sampleID), file("*trimming_report.txt"), emit: trim_stats

  script:
  log.info "----- Trim Galore Running on: ${sampleID} -----"

  if (params.non_directional) {
    directionality = '--non_directional'
  } 

  if (params.read_type == "SE"){
    paired_end = ''
  }

  if (params.read_type == "PE"){
    paired_end = '--paired'
  }

  if (params.workflow == "rrbs"){
    rrbs_flag = '--rrbs'
  }

  """
    trim_galore --cores ${task.cpus} ${paired_end} ${rrbs_flag} ${directionality} --gzip --length ${params.trimLength} -q ${params.qualThreshold}  --stringency ${params.adapOverlap}  -a ${params.adaptorSeq}  --fastqc ${fq_reads}
  """
}