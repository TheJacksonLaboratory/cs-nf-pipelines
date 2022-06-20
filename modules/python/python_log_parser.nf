process LOG_PARSER {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'logparser' }", pattern: "*.summary_QC_metrics.txt", mode: 'copy'

  container 'docker://python:3.8.10'

  input:
  tuple val(sampleID), file(log_cutadapt)
                       file(log_bowtie)
                       file(log_sorted_metrics)
                       file(log_mtdna_content)
                       file(log_pbc_qc)
                       file(log_fraction_reads)

  output:
  tuple val(sampleID), file("*.summary_QC_metrics.txt")

  script:
  log.info "----- LogParser on ${sampleID} -----"
  """
  python ${projectDir}/bin/atac/LogParser.py > ${sampleID}.summary_QC_metrics.txt
  """
}
