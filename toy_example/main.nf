#!/usr/bin/env nextflow

// 1.) enable DSL2
nextflow.enable.dsl=2

log.info """\
 TOY EXAMPLE   P I P E L I N E
 ===================================
 fq_path        : ${params.fq_path}
 outdir         : ${params.outdir}
 """

read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}${params.extension}",checkExists:true )

process trim {
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

// 5. Use RSEM for quantification 

process rsem_ref_pull {
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

process rsem_ref_build {
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

process rsem_expression {
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

workflow{
  ref_files=rsem_ref_build(rsem_ref_pull())
  trim_ch=trim(read_ch)
  rsem_expression(trim_ch, ref_files)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Output in : $params.outdir\n" : "Oops .. something went wrong" )
}
