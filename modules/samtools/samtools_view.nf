process SAMTOOLS_VIEW {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'samtools_view' }", pattern:"*.bam", mode:'copy', enabled: params.keep_intermediate

  input:
      tuple val(sampleID), file(sam)
      val(view_string)
      val(filename)

  output:
      tuple val(sampleID), file("*.bam"), emit: bam

  script:

    output_name = "${filename}" == 'emase' ? "${sam.baseName}.bam" : "${sampleID}_${filename}.bam"

    """
    samtools view ${view_string} ${sam} >  ${output_name}
    """
  
  stub:
    """
    ${sam.baseName}.bam 
    """
}
