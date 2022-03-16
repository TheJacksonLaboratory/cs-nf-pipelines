process RSEM_ALIGNMENT_EXPRESSION {
  tag "$sampleID"

  cpus 12
  memory { 60.GB * task.attempt }
  time { 24.h * task.attempt }
  clusterOptions '-q batch'
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

  if (params.read_prep == "stranded"){
    prob="--forward-prob 0"
  }
  if (params.read_prep == "non_stranded"){
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
// Toy Example RSEM below
process RSEM_REF_PULL {
  publishDir "${params.pubdir}/rsem/ref"

  output:
  tuple file("*.gtf"), file("*.fa")

  when:
  params.ref_pull=='true'

  script:
  log.info "----- RSEM Pull Genome Running on: ${sampleID} -----"
  """
  wget ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
  wget ftp://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/Mus_musculus.GRCm38.82.chr.gtf.gz
  gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz
  gunzip Mus_musculus.GRCm38.82.chr.gtf.gz
  """
}

process RSEM_REF_BUILD {
  publishDir "${params.pubdir}/rsem/ref"
  container "quay.io/jaxcompsci/rsem_bowtie2_star:0.1.0"

  input:
  tuple file(gtf), file(fa)

  output:
  file("*")

  script:
  log.info "----- RSEM Build Reference Running on: ${sampleID} -----"
        """
        rsem-prepare-reference \
        --gtf ${gtf} \
        --bowtie2 \
        ${fa} \
        ${params.species}

        """
}

process RSEM_EXPRESSION {
  publishDir "${params.pubdir}/rsem/exp"
  container "quay.io/jaxcompsci/rsem_bowtie2_star:0.1.0"

  input:
  tuple val(sampleId), file(R1), file(R2)
  file(ref_files)

  output:

  file "*"

  script:
  log.info "----- RSEM Expression Running -----"
  """
  rsem-calculate-expression -p 8 --paired-end \
  --bowtie2 \
  --estimate-rspd \
  --append-names \
  --output-genome-bam \
  ${R1} ${R2} \
  ${params.species} \
  Toy_Ex
  """
}

process RSEM_SIMULATE_READS{
  publishDir "${params.pubdir}/rsem/sim"
  container "quay.io/jaxcompsci/rsem_bowtie2_star:0.1.0"

  input:
  tuple file(estimated_model_file), file(estimated_isoform_results)

  output:
  file "*"

  script:
  log.info "----- RSEM Simulate Reads Running -----"
  """
  reference_name ${estimated_model_file} ${estimated_isoform_results} 0.2 50000000 simulated_reads
  """
}
