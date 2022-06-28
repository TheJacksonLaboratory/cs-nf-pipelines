process RSEM_ALIGNMENT_EXPRESSION {
  tag "$sampleID"

  cpus 12
  memory { 60.GB * task.attempt }
  time { 24.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

    container 'quay.io/jaxcompsci/rsem_bowtie2_star:0.1.0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'rsem' }", pattern: "*stats", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'rsem' }", pattern: "*results*", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'rsem' }", pattern: "*genome.bam", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'rsem' }", pattern: "*transcript.bam", mode:'copy'

  input:
  tuple val(sampleID), file(reads)
  file(rsem_ref_files)

  output:
  file "*stats"
  file "*results*"
  tuple val(sampleID), file("rsem_aln_*.stats"), emit: rsem_stats
  tuple val(sampleID), file("*genes.results"), emit: rsem_genes
  tuple val(sampleID), file("*isoforms.results"), emit: rsem_isoforms
  tuple val(sampleID), file("*.genome.bam"), emit: bam
  tuple val(sampleID), file("*.transcript.bam"), emit: transcript_bam

  script:
  log.info "----- Genome Alignment Running on: ${sampleID} -----"

  if (params.read_prep == "reverse_stranded") {
    prob="--forward-prob 0"
  }

  if (params.read_prep == "forward_stranded") {
    prob="--forward-prob 1"
  }

  if (params.read_prep == "non_stranded") {
    prob="--forward-prob 0.5"
  }

  if (params.read_type == "PE"){
    frag=""
    stype="--paired-end"
    trimmedfq="${reads[0]} ${reads[1]}"
  }
  if (params.read_type == "SE"){
    frag="--fragment-length-mean 280 --fragment-length-sd 50"
    stype=""
    trimmedfq="${reads[0]}"
  }
  """
  rsem-calculate-expression -p $task.cpus \
  ${prob} \
  ${stype} \
  ${frag} \
  --${params.rsem_aligner} \
  --append-names \
  --seed-length ${params.seed_length} \
  --output-genome-bam \
  ${trimmedfq} \
  ${params.rsem_ref_prefix} \
  ${sampleID} \
  2> rsem_aln_${sampleID}.stats
  """
}