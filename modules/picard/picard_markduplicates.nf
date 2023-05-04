process PICARD_MARKDUPLICATES {
  tag "$sampleID"

  cpus 1
  memory 16.GB
  time '12:00:00'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  // save if mouse and wes or save if keep intermediate
  publishDir {
      def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : '' 
      "${params.pubdir}/${ params.organize_by=='sample' ? type+sampleID+'/bam' : 'picard'}"
  }, pattern: "*.bam", mode: 'copy', enabled: params.gen_org=='mouse' ? true : params.keep_intermediate

  publishDir {
      def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
      "${params.pubdir}/${ params.organize_by=='sample' ? type+sampleID+'/stats' : 'picard'}"
  }, pattern: "*.txt", mode: 'copy'


  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*_dedup.bam"), emit: dedup_bam
  tuple val(sampleID), file("*_dedup.bai"), emit: dedup_bai
  tuple val(sampleID), file("*.txt"), emit: dedup_metrics

  script:
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  prefix = "${sampleID}.mLb.mkD"

  if (params.workflow == "atac"){
    output = "${sampleID}.sorted.marked4_dedup.bam"
  }
  if (params.workflow == "chipseq"){
    output = "${prefix}.sorted.marked4_dedup.bam"
  } 

  if (params.workflow != "atac" && params.workflow != "chipseq")
  """
  picard -Xmx${my_mem}G MarkDuplicates \
  I=${bam} \
  O=${sampleID}_dedup.bam \
  M=${sampleID}_dup_metrics.txt \
  REMOVE_DUPLICATES=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT
  """
  else
  """
  picard -Xmx${my_mem}G MarkDuplicates \
  I=${bam[0]} \
  O=${output} \
  M=${sampleID}.sorted.metrics.txt \
  REMOVE_DUPLICATES=false \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  TMP_DIR=${params.tmpdir} \
  > ${sampleID}.picard.log 2>&1  
  """
}
