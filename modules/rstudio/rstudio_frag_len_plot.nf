process FRAG_LEN_PLOT {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'rstudio' }", pattern: "*fraglen_plot.pdf", mode: 'copy'
  container 'quay.io/jaxcompsci/rstudio:4.2.0' 

  input:
  tuple val(sampleID), file(frag_len_count)

  output:
  tuple val(sampleID), file("*fraglen_plot.pdf")

  script:
  log.info "----- Fragment Length Plot on ${sampleID} -----"
  """
  Rscript ${projectDir}/bin/atac/fragment_length_plot.R ${frag_len_count}
  mv fraglen_plot.pdf ${sampleID}_fraglen_plot.pdf
  """
}
