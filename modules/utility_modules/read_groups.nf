process READ_GROUPS {
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'read_groups' }", pattern: "*read_group.txt", mode:'copy', enabled: params.workflow == 'rnaseq' || params.keep_intermediate

  input:
  tuple val(sampleID), file(fq_reads)
  val(picard)

  output:
  tuple val(sampleID), file("*.txt"), emit: read_groups

  script:
  if (picard=="picard"){
    p='-p'
  }
  else{
    p=''
  }
  """
  python ${projectDir}/bin/shared/read_group_from_fastq.py $p -s ${sampleID} -o ${sampleID}_read_group.txt ${fq_reads[0]}
  """
  }
