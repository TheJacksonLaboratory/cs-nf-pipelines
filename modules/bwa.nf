process BWA_MEM {
  tag "sampleID"

  cpus 1
  memory 5.GB
  time '02:00:00'
  clusterOptions '-q batch'

  container 'bwa-0.7.9a_python_2.7.3.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bwa_mem' }", pattern: "*.sam", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)
  tuple val(sampleID), file(read_groups)


  output:
  tuple val(sampleID), file("*.sam"), emit: bwa_mem

  script:
  log.info "----- BWA-MEM Alignment Running on: ${sampleID} -----"

  if (params.reads == "SE"){
    inputfq="${fq_reads[0]}"
    }
  if (params.reads == "PE"){
    inputfq="${fq_reads[0]} ${fq_reads[1]}"
    }

  """
  rg=\$(cat ${read_groups})
  bwa mem -M -B ${params.mismatch_penalty} -t ${task.cpus} -R \${rg} ${params.ref_fa} $inputfq > ${sampleID}.sam
  """
  }
