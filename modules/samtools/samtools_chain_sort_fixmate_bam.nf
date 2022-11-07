process CHAIN_SORT_FIXMATE_BAM {
  tag "$sampleID"

  cpus  8
  memory 20.GB
  time '20:00:00'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'samtools' }", pattern: "*.filtered.shifted.*", mode: 'copy'
  container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

  input:
  tuple val(sampleID), file(bam_mm10)

  output:
  tuple val(sampleID), file("*.filtered.shifted.*")

  when: params.chain != null

  script:
  // This module is for Non-Reference Strain Samples. 
  // To sort bam by read name, fix the mate information, re-sort by coordinates and filter Mitochondrial Reads from bam file. 
  """
  # sort bam by read name
  samtools sort \
  -n \
  -@ $task.cpus -O bam \
  -o ${sampleID}.tmp3.mm10.bam ${bam_mm10[0]}

  # fix the mate information. This is done to fix 'TLEN' which is required for MACS2
  samtools fixmate \
  -O bam ${sampleID}.tmp3.mm10.bam ${sampleID}.tmp4.mm10.bam

  # re-sort bam by coordinates
  samtools sort \
  -@ $task.cpus -O bam \
  -o ${sampleID}.tmp5.mm10.bam ${sampleID}.tmp4.mm10.bam

  samtools index ${sampleID}.tmp5.mm10.bam

  # filter Mitochondrial Reads
  samtools view ${sampleID}.tmp5.mm10.bam \
  -h 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
  | grep -ve 'SN:MT*' > tmp.sam

  samtools view -b tmp.sam \
  > ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.shifted.mm10.bam

  samtools index \
  ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.shifted.mm10.bam
  """
}
