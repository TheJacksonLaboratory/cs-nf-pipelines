
process TRIM {
  publishDir "${params.outdir}/trimmed"

  input:
  tuple val(sampleId), file(reads)

  output:
  tuple val(sampleId), file('*R1_paired.fastq.gz'), file('*R2_paired.fastq.gz')

  script:
  """
  trimmomatic \
  PE \
  ${params.fq_path}/${reads[0]} \
  ${params.fq_path}/${reads[1]} \
  ${sampleId}_R1_paired.fastq.gz \
  ${sampleId}_R1_unpaired.fastq.gz \
  ${sampleId}_R2_paired.fastq.gz \
  ${sampleId}_R2_unpaired.fastq.gz \
  LEADING:${params.t_lead} \
  TRAILING:${params.t_trail} \
  MINLEN:${params.min_len}
  """

}

// 3. Use RSEM for quantification:
// Three parts: pull reference genome, make reference files, then quantify

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

  input:
  tuple file(gtf), file(fa)

  when:
  params.ref_build=='true'

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
