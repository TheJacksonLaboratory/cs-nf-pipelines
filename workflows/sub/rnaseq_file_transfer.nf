process TRANSFER_FILES_HSA {

  tag "sampleID"

  input:
  tuple sampleID, file(rsem_stat)
  tuple sampleID, file(sum_stats)
  tuple sampleID, file(fastq_stat)
  tuple sampleID, file(gene_results)
  tuple sampleID, file(picard_metrics)
  tuple sampleID, file(isoform_results)
  tuple sampleID, file(avg_med_coverage)

  when:
  params.gen_org == "human"

  script:
  log.info "----- Moving Files to Output Directory for: ${sampleID} -----"

  """
  echo "${sample_tmpdir}" > sampletmpdir.txt
  fullpath=`cat sampletmpdir.txt`
  finfullpath=\$(basename \$fullpath)
  mv "\${fullpath}_tmp" "${params.outdir}/\${finfullpath}/"
  """
  }

process TRANSFER_FILES_MMU {

  tag "sampleID"

  input:
  tuple sampleID, file(rsem_stat)
  tuple sampleID, file(fastq_stat)
  tuple sampleID, file(gene_results)
  tuple sampleID, file(picard_metrics)
  tuple sampleID, file(isoform_results)


  when:
  params.gen_org == "mouse"

  script:
  log.info "----- Moving Files to Output Directory for: ${sampleID} -----"

  """
  echo "${sample_tmpdir}" > sampletmpdir.txt
  fullpath=`cat sampletmpdir.txt`
  finfullpath=\$(basename \$fullpath)
  mv "\${fullpath}_tmp" "${params.outdir}/\${finfullpath}/"
  """
}
