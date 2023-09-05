process FASTP {

  tag "$sampleID"

  cpus 6
  memory 64.GB
  time '24:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
  
  container 'quay.io/biocontainers/fastp:0.23.2--h5f740d0_3'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'stats'}", pattern: "${sampleID}_fastp_report.*", mode:'copy'

  input:
    tuple val(sampleID), file(fq_reads)

  output:
    tuple val(sampleID), file("${sampleID}_fastp_report.html"), emit: quality_html
    tuple val(sampleID), file("${sampleID}_fastp_report.json"), emit: quality_json
    tuple val(sampleID), file("${sampleID}.trimmed.R*.fastq"), emit: trimmed_fastq

  script:

  if (!params.fastq2)
    """
    fastp -i ${fq_reads[0]} \
          -o ${sampleID}.trimmed.R1.fastq \
          -q ${params.quality_phred} \
          -u ${params.unqualified_perc} \
          -w ${task.cpus} \
          -j ${sampleID}_fastp_report.json \
          -h ${sampleID}_fastp_report.html \
          -R "${sampleID} fastp report"
    """
  else
    """
    fastp -i ${fq_reads[0]} \
          -I ${fq_reads[1]} \
          -o ${sampleID}.trimmed.R1.fastq \
          -O ${sampleID}.trimmed.R2.fastq \
          -q ${params.quality_phred} \
          -u ${params.unqualified_perc} \
          -w ${task.cpus} \
          -j ${sampleID}_fastp_report.json \
          -h ${sampleID}_fastp_report.html \
          -R "${sampleID} fastp report"
    """
}