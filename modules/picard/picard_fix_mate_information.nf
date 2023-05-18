process PICARD_FIX_MATE_INFORMATION {
  tag "$sampleID"

  cpus = 1
  memory { bam.size() < 30.GB ? 6.GB : 48.GB }
  time { bam.size() < 30.GB ? '03:00:00' : '24:00:00' }
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
  
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*fixed_mate.bam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*_fixed_mate.bam"), emit: fixed_mate_bam

  script:

  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  picard -Xmx${my_mem}G FixMateInformation \
  I=${bam} \
  O=${sampleID}_fixed_mate.bam \
  TMP_DIR=${workDir}/temp \
  ADD_MATE_CIGAR=true
  """
}
