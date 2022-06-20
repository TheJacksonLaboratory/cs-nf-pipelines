process FLAG_PCR_DUPES {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.sorted.metrics", mode: 'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.log", mode: 'copy'
  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(sampleID), file(bam_file)
  tuple val(sampleID), file(bai_file)

  output:
  tuple val(sampleID), file("*.sorted.marked.bam"), emit: marked_bam
  tuple val(sampleID), file("*.sorted.marked.bai"), emit: marked_bai
  tuple val(sampleID), file("*.sorted.metrics"), emit: srt_metr_log
  tuple val(sampleID), file("*.picard.log")

  script:
  log.info "----- Flagging PCR Duplicates on : ${sampleID} -----"
  """
  gatk \
  --java-options "-XX:ParallelGCThreads=$task.cpus -Djava.io.tmpdir=${params.tmpdir}" \
  MarkDuplicates \
  --INPUT ${bam_file} \
  --OUTPUT ${sampleID}.sorted.marked.bam \
  --METRICS_FILE ${sampleID}.sorted.metrics \
  --REMOVE_DUPLICATES false \
  --CREATE_INDEX true \
  --VALIDATION_STRINGENCY LENIENT \
  --TMP_DIR ${params.tmpdir} \
  > ${sampleID}.picard.log 2>&1
  """

}
