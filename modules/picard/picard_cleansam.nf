process PICARD_CLEANSAM {
  tag "$sampleID"

  cpus = 1
  memory { bam.size() < 60.GB ? 8.GB : 24.GB }
  time { bam.size() < 60.GB ? '06:00:00' : '12:00:00' }

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*_cleaned.bam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*_cleaned.bam"), emit: cleaned_bam

  script:

  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  picard -Xmx${my_mem}G CleanSam \
  I=${bam} \
  TMP_DIR=${workDir}/temp \
  O=${sampleID}_cleaned.bam
  """
}
