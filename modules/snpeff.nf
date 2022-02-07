process SNPEFF{
  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  // try with old container
  container 'gatk-3.6_snpeff-3.6c_samtools-1.3.1_bcftools-1.11.sif'
// this is most recent but does not accept the old .bin files (v4 did not work either)
//  container 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'snpeff' }", pattern:"*.*", mode:'copy'

  input:
  tuple val(sampleID),file(vcf)

  output:
  tuple val(sampleID),file("*.vcf"), emit:vcf
  tuple val(sampleID),file("*.html")

  // may change -v to a paramiter
  script:
  log.info "----- snpEff Running on: ${sampleID} -----"
  """
  java -Xmx8g -jar /snpEff/snpEff.jar GRCm38.75 \
  -c ${params.snpEff_config} \
  -o gatk \
  -s ${sampleID}_snpeff.html \
  ${vcf} > ${sampleID}_snpeff.vcf
  """
}
process SNPEFF_HUMAN{
  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  // try with old container
  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'
// this is most recent but does not accept the old .bin files (v4 did not work either)
//  container 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'snpeff' }", pattern:"*.*", mode:'copy'

  input:
  tuple val(sampleID),file(vcf)
  val(indel_snp)

  output:
  tuple val(sampleID),file("*.vcf"), emit:vcf

  script:
  log.info "----- snpEff Running on: ${sampleID} -----"
  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /snpEff_v4_3/snpEff/snpEff.jar \
  -v -lof ${params.gen_ver} \
  -dataDir ${params.hgvs_data} \
  -noStats ${vcf}  > ${sampleID}_snpeff_${indel_snp}.vcf
  """
}
