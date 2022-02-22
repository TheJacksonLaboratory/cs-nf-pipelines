process MAKE_VCF_LIST {
 
  input:
  tuple val(sampleID), val(chromes)

  output:
  tuple val(sampleID), file("*.list"), emit: list

  script:
  log.info "----- Make VCF List from Chromosomes: ${sampleID} -------"

  // Puts Individual Chromosome Files into List for MergeVCFs 
  def my_list=chromes.join("\n")

  """
  echo "$my_list" > ${sampleID}.list
  """
}
