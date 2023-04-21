process QUALITY_CHECKS {
  tag "$sampleID"

  cpus 2
  memory 4.GB
  time '04:00:00'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats'  : 'samtools' }", pattern: "*.fragment_length_count.txt", mode: 'copy'
  container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

  input:
  tuple val(sampleID), file(sort_rm_filter_bam)

  output:
  tuple val(sampleID), file("*.fragment_length_count.txt")

  script:
  // Get the fragment length count from bam file for Quality Checks.
  """
  samtools view \
  -@ $task.cpus ${sort_rm_filter_bam[0]} \
  | awk '\$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n \
  | sed -e 's/^[ \\t]*//' > ${sampleID}.fragment_length_count.txt
  """
}
