process RSEM_ALIGNMENT_EXPRESSION {
  tag "$sampleID"

  cpus 12
  memory { 60.GB * task.attempt }
  time { 24.h * task.attempt }
  errorStrategy 'finish'
  maxRetries 1

    container 'quay.io/jaxcompsci/rsem_bowtie2_star:0.1.0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'rsem' }", pattern: "*stats", mode:'copy', enabled: params.rsem_aligner == "bowtie2"
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'rsem' }", pattern: "*results*", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'rsem' }", pattern: "*genome.sorted.ba*", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'rsem' }", pattern: "*transcript.sorted.ba*", mode:'copy'

  input:
  tuple val(sampleID), path(reads)
  path(rsem_ref_files)
  val(rsem_ref_prefix)

  output:
  path "*stats"
  path "*results*"
  tuple val(sampleID), path("rsem_aln_*.stats"), emit: rsem_stats
  tuple val(sampleID), path("*.stat/*.cnt"), emit: rsem_cnt
  tuple val(sampleID), path("*genes.results"), emit: rsem_genes
  tuple val(sampleID), path("*isoforms.results"), emit: rsem_isoforms
  tuple val(sampleID), path("*.genome.bam"), emit: bam
  tuple val(sampleID), path("*.transcript.bam"), emit: transcript_bam
  tuple val(sampleID), path("*.genome.sorted.bam"), path("*.genome.sorted.bam.bai"), emit: sorted_genomic_bam
  tuple val(sampleID), path("*.transcript.sorted.bam"), path("*.transcript.sorted.bam.bai"), emit: sorted_transcript_bam
 
  script:

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
  if (params.rsem_aligner == "bowtie2"){
    outbam="--output-genome-bam --sort-bam-by-coordinate"
    seed_length="--seed-length ${params.seed_length}"
  }
  if (params.rsem_aligner == "star") {
    outbam="--star-output-genome-bam --sort-bam-by-coordinate"
    seed_length=""
  }

  """
  rsem-calculate-expression -p $task.cpus \
  ${prob} \
  ${stype} \
  ${frag} \
  --${params.rsem_aligner} \
  --append-names \
  ${seed_length} \
  ${outbam} \
  ${trimmedfq} \
  ${rsem_ref_prefix} \
  ${sampleID} \
  2> rsem_aln_${sampleID}.stats
  """
}