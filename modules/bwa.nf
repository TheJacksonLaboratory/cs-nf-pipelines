process BWA_MEM {
  tag "sampleID"

  cpus 4
  memory 16.GB
  time '30:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/bwa:0.7.3a--h5bf99c6_6'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bwa_mem' }", pattern: "*.sam", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)
  tuple val(sampleID), file(read_groups)

  output:
  tuple val(sampleID), file("*.sam"), emit: sam

  script:
  log.info "----- BWA-MEM Alignment Running on: ${sampleID} -----"

  if (params.read_type == "SE"){
    inputfq="${fq_reads[0]}"
    }
  if (params.read_type == "PE"){
    inputfq="${fq_reads[0]} ${fq_reads[1]}"
    }

  """
  rg=\$(cat $read_groups)
  rg=\$(echo \$rg | sed -r 's/[=]+/:/g')
  rg=\$(echo \$rg | sed -r 's/[RG]+//g')
  rg=\$(echo \$rg | sed -r 's/[ ]+/\t/g')

  echo '\'\$rg\''

  bwa mem -M -R "@RG\t\${rg}" \
  ${params.ref_fa} $inputfq > ${sampleID}.sam
  """
  }
