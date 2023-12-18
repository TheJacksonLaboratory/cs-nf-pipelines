process BWA_MEM_HLA {
  tag "$sampleID"

  cpus 8
  memory 60.GB
  time 48.hour
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bwa_mem' }", pattern: "*.bam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), path(fq_reads), file(read_groups)

  output:
  tuple val(sampleID), path("*.bam"), emit: bam

  script:

  if (params.read_type == "SE"){
    inputfq="${fq_reads[0]}"
    }
  if (params.read_type == "PE"){
    inputfq="${fq_reads[0]} ${fq_reads[1]}"
    }

  """
  rg=\$(cat $read_groups)

  run-bwamem -t $task.cpus -R \${rg} -o ${sampleID} -H ${params.ref_fa_indices} $inputfq | sh

  """
}