process BAMTOOLS_RNASEQ_MOUSE {

  tag "sampleID"

  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/bamtools:2.5.1--h9a82719_9'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'quality_stats' }", pattern:"*.txt", mode:'copy'

  input:
  tuple val(sampleID), file(reordered_sorted_bam)

  output:
  tuple val(sampleID), file("*metrics.txt"), emit: picard_metrics

  script:
  if (params.read_type == "PE")

    """
    bamtools stats -insert -in ${reordered_sorted_bam} > ${sampleID}_aln_metrics.txt
    """

  else if (params.read_type == "SE")

    """
    bamtools stats -in ${reordered_sorted_bam} > ${sampleID}_aln_metrics.txt
    """
}
