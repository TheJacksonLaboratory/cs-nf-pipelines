// RNASEQ CAROLYN
process RSEM_ALIGNMENT_EXPRESSION {

  tag "sampleID"

  cpus 12
  memory 30.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'rsem_bowtie2_samtools_picard.v2.sif'

  if (${params.organize_by}=='analysis'){
    publishDir "${params.pubdir}/rsem", pattern: "*stats", mode: 'copy'
    publishDir "${params.pubdir}/rsem", pattern: "*results*", mode: 'copy'
    publishDir "${params.pubdir}/rsem", pattern: "*.bam", mode: 'copy'
  }
  else if (${params.organize_by}=='sample'){
    publishDir "${params.pubdir}/${sampleID}", pattern: "*stats", mode: 'copy'
    publishDir "${params.pubdir}/${sampleID}", pattern: "*results*", mode: 'copy'
    publishDir "${params.pubdir}/${sampleID}", pattern: "*.bam", mode: 'copy'
  }

  input:
  tuple val(sampleID), file(reads)
  file(rsem_ref_files)

  output:
  file "*stats"
  file "*results*"
  tuple val(sampleID), file("rsem_aln_*.stats"), emit: rsem_stats
  tuple val(sampleID), file("*genes.results"), emit: rsem_genes
  tuple val(sampleID), file("*isoforms.results"), emit: rsem_isoforms
  tuple val(sampleID), file("*genome.bam"), emit: genome_sorted_bam


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
  rsem-calculate-expression -p 12 \
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
  """
  wget ftp://ftp.ensembl.org/pub/release-82/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
  wget ftp://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/Mus_musculus.GRCm38.82.chr.gtf.gz
  gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz
  gunzip Mus_musculus.GRCm38.82.chr.gtf.gz
  """
}

process RSEM_REF_BUILD {
  publishDir "${params.pubdir}/rsem/ref"
  container "dceoy/rsem"

  input:
  tuple file(gtf), file(fa)

  output:
  file("*")

  script:
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
  container "dceoy/rsem"

  input:
  tuple val(sampleId), file(R1), file(R2)
  file(ref_files)

  output:
  file "*"

  script:
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
  container "dceoy/rsem"

  input:
  tuple file(estimated_model_file), file(estimated_isoform_results)

  output:
  file "*"

  script:
  """
  reference_name ${estimated_model_file} ${estimated_isoform_results} 0.2 50000000 simulated_reads
  """
}
