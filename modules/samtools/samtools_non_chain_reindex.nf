process NON_CHAIN_REINDEX {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'samtools' }", pattern: "*.${params.processed_bam_extension}.ba*", mode: 'copy'
  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(bam_shifted)

  output:
  tuple val(sampleID), file("*.${params.processed_bam_extension}.ba*")

  when: params.chain == null

  script:
  log.info "----- Performing name sort the bam on ${sampleID} -----"
  """
  samtools index ${bam_shifted[0]}

  samtools view ${bam_shifted[0]} \
  -h 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
  | grep -ve 'SN:MT*\\|SN:GL*\\|SN:JH:*' > tmp.sam

  samtools view -b tmp.sam \
  > ${sampleID}.${params.processed_bam_extension}.bam

  samtools index \
  ${sampleID}.${params.processed_bam_extension}.bam
  """
}
