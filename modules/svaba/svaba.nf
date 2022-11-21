process SvABA {
  tag "$meta.patient"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'

  container 'quay.io/jaxcompsci/svaba::latest'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'svaba' }", pattern: "*.log", mode:'copy'

  input:
  tuple val(sampleID), file(normal_bam), file(normal_bai), val(meta)
  tuple val(sampleID), file(tumor_bam), file(normal_bai), val(meta)

  output:
  tuple val(meta), file("*_svaba.vcf"), emit: svaba_vcf
  tuple val(meta), file("*_svaba.bps.txt.gz"), emit: svaba_unfiltered_variants
  tuple val(meta), file("*_svaba.contigs.bam"), emit: svaba_contigs_bam
  tuple val(meta), file("*_svaba.discordants.txt.gz"), emit: svaba_discordants
  tuple val(meta), file("*_svaba.log"), emit: svaba_log
  tuple val(meta), file("*_svaba.alignments.txt.gz"), emit: svaba_alignments

  script:
  log.info "----- SvABA Running on: ${sampleID} -----" 
  //Estimate somatic variants using SvABA
  
  """
  svaba run -t ${tumor_bam} \
    -n ${normal_bam} \
    -p ${cpus} \
    -a ${meta.patient}_svaba \
    -G ${params.ref_fa} 
  """