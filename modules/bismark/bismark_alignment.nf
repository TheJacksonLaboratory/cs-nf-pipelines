process BISMARK_ALIGNMENT {
  tag "$sampleID"

  cpus 8
  memory {60.GB * task.attempt}
  time {30.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'quay.io/biocontainers/bismark:0.23.1--hdfd78af_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/alignment' : 'bismark_align' }", pattern: "*.bam", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'bismark_align' }", pattern: "*report.txt", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/alignment' : 'bismark_align' }", pattern: "*unmapped*", mode:'copy', enabled: params.keep_intermediate
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/alignment' : 'bismark_align' }", pattern: "*ambiguous*", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*.bam"), emit: bam
  tuple val(sampleID), file("*report.txt"), emit: report
  tuple val(sampleID), file("*ambiguous*"), emit: ambiguous_reads
  tuple val(sampleID), file("*unmapped*"), emit: unmapped_reads

  script:
  log.info "----- Bismark Alignment Running on: ${sampleID} -----"

  inputfq = params.read_type == 'PE' ?  "-1 ${fq_reads[0]} -2 ${fq_reads[1]}" : "-1 ${fq_reads[0]}"
  directionality = params.non_directional ? '--non_directional': ''

  aligner = params.aligner == "bismark_hisat" ? "--hisat2" : "--bowtie2"

  """
  bismark ${aligner} --bam -p ${task.cpus} ${directionality} -L ${params.seedLength} -N ${params.seedMismatch} -minins ${params.MinInsert} -maxins ${params.MaxInsert}  --unmapped --ambiguous --genome ${params.ref_fa_index} ${inputfq}
  """
}
