process SHORT_ALIGNMENT_MARKING {
  tag "$sampleID"

  cpus 1
  memory 24.GB
  time '24:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'short_alignment_marking' }", pattern:"*.marked.bam", mode:'copy', enabled: params.keep_intermediate
  
  input:
  tuple val(sampleID), file(aligned_bam)

  output:
  tuple val(sampleID), file("*.marked.bam"), emit: marked_bam

  script:
  // parses the bam file and marks as unmapped a read with alignment length below a user-defined threshold. Reads are not filtered from the bam file but kept as unmapped.
  """
  ${projectDir}/bin/pta/filter_bam -I ${aligned_bam} -A1 30 -A2 30 -o ${sampleID}.marked.bam | samtools view -b -o ${sampleID}.marked.bam
  """
}

/*
-A1, --ALN_LEN_PRIM ALN_LEN_PRIM
     Primary (loose) alignment length
-A2, --ALN_LEN_SECOND ALN_LEN_SECOND
     Supplementary (strict) alignment length
-o, --OUT_PREFIX OUT_PREFIX
    Output file prefix

NOTE: -o does not actually produce an output file. 
NOTE: The BAM file produced here, is corrupt. It requires sorting and cleaning (non mapped reads have non 0 MAPQ) and mate information to be fixed. 
*/