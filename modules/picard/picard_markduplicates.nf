process PICARD_MARKDUPLICATES {
  tag "$sampleID"

  cpus 1
  memory { bam.size() < 60.GB ? 16.GB : 32.GB }
  time { bam.size() < 60.GB ? '12:00:00' : '24:00:00' }

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  // save if mouse or save if keep intermediate
  publishDir {
      def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : '' 
      "${params.pubdir}/${ params.organize_by=='sample' ? type+sampleID+'/bam' : 'picard'}"
  }, pattern: "*.bam", mode: 'copy', enabled: params.gen_org=='mouse' || params.workflow=='chipseq' ? true : params.keep_intermediate

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

  """
  picard -Xmx${my_mem}G MarkDuplicates \
  I=${bam[0]} \
  O=${sampleID}_dedup.bam \
  M=${sampleID}_dup_metrics.txt \
  REMOVE_DUPLICATES=false \
  CREATE_INDEX=true \
  TMP_DIR=${workDir}/temp \
  VALIDATION_STRINGENCY=SILENT
  """
}
