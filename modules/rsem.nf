// RNASEQ CAROLYN
process RSEM_ALIGNMENT_EXPRESSION {

  tag "sampleID"

  cpus 12
  memory 30.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'rsem_bowtie2_samtools_picard.v2.sif'

  // ? are these cleared out? how necessary?
  publishDir "${sample_tmpdir}_tmp", pattern: "*stats", mode: 'copy'
  publishDir "${sample_tmpdir}_tmp", pattern: "*results*", mode: 'copy'

  input:
  tuple val(sampleID), file(trimmed)

  output:
  file "*stats"
  file "*results*"
  tuple sampleID, file("*genome.bam")
  tuple sampleID, file("*aln.stats")
  tuple sampleID, file("*genes.results")
  tuple sampleID, file("*isoforms.results")

  script:
  log.info "-----Genome Alignment Running on: ${sampleID} -----"

  if (params.read_prep == "stranded"){
    prob="--forward-prob 0"
  }
  if (params.read_prep == "non_stranded"){
    prob="--forward-prob 0.5"
  }

  if (params.reads == "PE"){
    frag=""
    trimmedfq="--paired-end ${trimmed[0]} ${trimmed[1]}"
  }
  if (params.reads == "SE"){
    frag="--fragment-length-mean 280 --fragment-length-sd 50"
    trimmedfq="${trimmed[0]}"
  }

  """
  rsem-calculate-expression -p 12 \
  --phred33-quals $frag \
  --seed-length ${params.seed_length} $prob \
  --time \
  --output-genome-bam ${params.aligner} \
  $trimmedfq ${params.rsem_ref_prefix} ${sampleID} \
  2> ${sampleID}_rsem_aln.stats
  """
}

process RSEM_REF_PULL {
  publishDir "${params.outdir}/rsem/ref"

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
  publishDir "${params.outdir}/rsem/ref"
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
  publishDir "${params.outdir}/rsem/exp"
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
  publishDir "${params.outdir}/rsem/sim"
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
