process CHAIN_BAD2UNIQ_READS {
  tag "$sampleID"

  cpus = 1

  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(bad_reads)

  output:
  tuple val(sampleID), file("ReadName_unique"), emit: uniq_reads
  
  when: params.chain != null

  shell:
  log.info "----- Removing 'bad reads' from bam file on ${sampleID} -----"
  '''
  cat !{bad_reads} \
  | awk '{print $5}' \
  | sed -r 's/\\,//g' \
  | sort -n \
  | uniq -c \
  | awk '{print $2}' \
  > ReadName_unique
  '''
}
