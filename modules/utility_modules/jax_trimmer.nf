process JAX_TRIMMER {

  tag "$sampleID"

  cpus 1
  memory 30.GB
  time '30:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'quality_stats' }", pattern: "*_stat", mode:'copy'

  input:
  tuple val(sampleID), path(fq_reads)

  output:
  tuple val(sampleID), file("*_stat"), emit: quality_stats
  tuple val(sampleID), file("*filtered_trimmed"), emit: trimmed_fastq

  script:

  if (params.read_type == "SE"){
    mode_HQ="-S -M"
    inputfq="${fq_reads[0]}"
  }
  if (params.read_type == "PE"){
    mode_HQ="-M"
    inputfq="${fq_reads[0]} ${fq_reads[1]}"
  }

  """
  python ${projectDir}/bin/shared/filter_trim.py $mode_HQ ${params.min_pct_hq_reads} -p ${params.hq_pct} $inputfq
  """
}
