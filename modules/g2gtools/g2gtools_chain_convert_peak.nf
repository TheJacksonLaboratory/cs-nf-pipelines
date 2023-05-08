process CHAIN_CONVERT {
  tag "$sampleID"

  cpus 1
  memory 10.GB
  time '10:00:00'

  container 'quay.io/jaxcompsci/g2gtools:0.1.31'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'g2gtools' }", pattern: "*.log", mode:'copy'

  input:
  tuple val(sampleID), file(bam_shifted)

  output:
  tuple val(sampleID), file("*.tmp.mm10.bam"), emit: coverted_bam
  tuple val(sampleID), file("*g2gconvert.log"), emit: log

  when: params.chain != null

  script:
  """
  g2gtools convert \
  -r -f bam -c ${params.chain} \
  -i ${bam_shifted[0]} \
  -o ${sampleID}.tmp.mm10.bam 2> ${sampleID}_g2gconvert.log
  """
}
