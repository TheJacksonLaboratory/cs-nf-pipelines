process SNPEFF{
  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  // SNPEFF and SNPSIFT need updating
  container 'gatk-3.6_snpeff-3.6c_samtools-1.3.1_bcftools-1.11.sif'
// this is most recent but does not accept the old .bin files (v4 did not work either)
// probably need to update snpEff downloadable files to update to newer version (v5.1)
//  container 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'snpeff' }", pattern:"*.*", mode:'copy'

  input:
  tuple val(sampleID),file(vcf)

  output:
  tuple val(sampleID),file("*.vcf"), emit:vcf
  tuple val(sampleID),file("*.html")
  // tuple val(sampleID),file("*")

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

  // SNPEFF and SNPSIFT need updating
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
process SNPEFF_ONEPERLINE {
  tag "sampleID"

  cpus 1
  memory 2.GB
  time '00:10:00'
  clusterOptions '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

  input:
  tuple val(sampleID), file(vcf)
  val(indel_snp)

  output:
  tuple	val(sampleID), file("*.vcf"), emit: vcf

  script:
  // the pl in here needs to be discovered. this will happen when making the container cook book
  """
  cat ${vcf} | /snpEff_v4_3/snpEff/scripts/vcfEffOnePerLine.pl > ${sampleID}_oneperline_${indel_snp}.vcf
  """
}
