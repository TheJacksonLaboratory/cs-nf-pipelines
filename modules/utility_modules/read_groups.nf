process READ_GROUPS {
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'

  container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'read_groups' }", pattern: "*read_group.txt", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(fq_reads)
  val(picard)

  output:
  tuple val(sampleID), file("*.txt"), emit: read_groups

  script:
  log.info "----- Read Group Information Determination Running on: ${sampleID} -----"
  if (picard=="picard"){
    p='-p'
  }
  else{
    p=''
  }
  """
  python ${projectDir}/bin/shared/read_group_from_fastq.py $p -o ${sampleID}_read_group.txt ${fq_reads[0]}
  """
  }
