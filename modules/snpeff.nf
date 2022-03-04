process SNPEFF{
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  // SNPEFF and SNPSIFT need updating
  container 'snpeff_5.1--hdfd78af_1.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'snpeff' }", pattern:"*.*", mode:'copy'

  input:
  tuple val(sampleID),file(vcf)

  output:
  tuple val(sampleID),file("*.vcf"), emit:vcf
  tuple val(sampleID),file("*.html")
  // tuple val(sampleID),file("*")

   script:
  log.info "----- snpEff Running on: ${sampleID} -----"
  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /usr/local/share/snpeff-5.1-1/snpEff.jar \
  ${params.gen_ver} \
  -c ${params.snpEff_config} \
  -o gatk \
  -s ${sampleID}_snpeff.html \
  -dataDir ${params.mvs_data} \
  ${vcf} > ${sampleID}_snpeff.vcf
  """
}
process SNPEFF_HUMAN{
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  container 'snpeff_5.1--hdfd78af_1.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'snpeff' }", pattern:"*.*", mode:'copy'

  input:
  tuple val(sampleID),file(vcf)
  val(indel_snp)

  output:
  tuple val(sampleID),file("*.vcf"), emit:vcf

  script:
  log.info "----- snpEff Running on: ${sampleID} -----"
  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /usr/local/share/snpeff-5.1-1/snpEff.jar  \
  -v -lof ${params.gen_ver} \
  -c ${params.snpEff_config} \
  -dataDir ${params.hgvs_data} \
  -noStats ${vcf}  > ${sampleID}_snpeff_${indel_snp}.vcf
  """
}
process SNPEFF_ONEPERLINE {
  tag "$sampleID"

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
