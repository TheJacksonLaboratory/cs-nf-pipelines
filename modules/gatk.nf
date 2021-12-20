// part A
process GATK_STATS_A {

  tag "sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'bedtools_2.27.1_python_2.7.3_java_1.8_GATK_3.4_samtools_1.3.1.sif'

  input:
  tuple sampleID, file(reord_sorted_bam)

  output:
  tuple sampleID, file("*gatk_temp3*"), emit: gatk_3
  tuple sampleID, file("*gatk_temp6*"), emit: gatk_6

  when:
  params.gen_org == "human"

  script:
  log.info "----- Human GATK Coverage Stats, Part 1 Running on: ${sampleID} -----"

  """
  samtools index ${reord_sorted_bam}

  java -Djava.io.tmpdir=$TMPDIR -Xmx48g -jar /GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${params.ref_fa} \
  --outputFormat table \
  -o ${sampleID}_gatk_temp1.txt \
  -I ${reord_sorted_bam} \
  -L  ${params.probes} \
  --omitPerSampleStats \
  --omitIntervalStatistics \
  --omitLocusTable \
  -U ALLOW_N_CIGAR_READS

  java -Djava.io.tmpdir=$TMPDIR -Xmx48g -jar /GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${params.ref_fa} \
  --outputFormat table \
  -o ${sampleID}_gatk_temp4.txt \
  -I ${reord_sorted_bam} \
  -L ${params.ctp_genes} \
  --omitPerSampleStats \
  --omitIntervalStatistics \
  --omitLocusTable \
  -U ALLOW_N_CIGAR_READS

  ${params.gatk_form} ${sampleID}_gatk_temp1.txt ${sampleID}_gatk_temp2.txt ${sampleID}_gatk_temp3.txt ${params.probes}

  ${params.gatk_form} ${sampleID}_gatk_temp4.txt ${sampleID}_gatk_temp5.txt ${sampleID}_gatk_temp6.txt ${params.ctp_genes}

  """
  }

// part B
process GATK_STATS_B {

  tag "sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'python_2_7_3'

  publishDir "${outdir}/gatk", pattern: "*.*", mode: 'copy'

  input:
  tuple sampleID, file(gatk_3)
  tuple sampleID, file(gatk_6)

  output:
  file "*CCP_interval_avg_median_coverage.bed"
  file "*exome_interval_avg_median_coverage.bed"
  tuple sampleID, file("*CP_interval_avg_median_coverage.bed")

  when:
  params.gen_org == "human"

  script:
  log.info "-----Human GATK Coverage Stats, Part 2 Running on: ${sampleID} -----"

  """

  python ${params.cov_calc} ${sampleID}_gatk_temp3.txt ${sampleID}_exome_interval_avg_median_coverage.bed

  python ${params.cov_calc} ${sampleID}_gatk_temp6.txt ${sampleID}_CCP_interval_avg_median_coverage.bed

  """
  }
