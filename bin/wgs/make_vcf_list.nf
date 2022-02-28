process MAKE_VCF_LIST {

  input:
  tuple val(sampleID), val(chroms)

  output:
  tuple val(sampleID), file("*.list"), emit: list

  script:
  log.info "----- Make VCF List from Chromosomes: ${sampleID} -------"

  // Puts Individual Chromosome Files into List for MergeVCFs
  def my_list=chroms.join("\n")

  // sort list by chromosome 1-x

  """
  echo "$my_list" > ${sampleID}.list
  """
}
