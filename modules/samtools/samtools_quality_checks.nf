process QUALITY_CHECKS {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'samtools' }", pattern: "*.fragment_length_count.txt", mode: 'copy'
  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(sort_rm_filter_bam)

  output:
  tuple val(sampleID), file("*.fragment_length_count.txt")

  script:
  log.info "----- Quality checks on ${sampleID} -----"
  log.info "----- Fragment/Insert size on ${sampleID} -----"
  """
  samtools view \
  -@ $task.cpus ${sort_rm_filter_bam[0]} \
  | awk '\$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n \
  | sed -e 's/^[ \\t]*//' > ${sampleID}.fragment_length_count.txt
  """
}
