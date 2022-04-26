process MAKE_VCF_LIST {
  tag "$sampleID"

  input:
  tuple val(sampleID), val(chroms)
  val(order)

  output:
  tuple val(sampleID), file("*.list"), emit: list

  script:
  log.info "----- Make VCF List from Chromosomes: ${sampleID} -------"

  // Puts Individual Chromosome Files In Order and Then Into List for MergeVCFs
  // convert paths to strings
  string_list = [] 
  for (int i = 0; i < chroms.size(); i++) {
    string_list.add(chroms[i].toString())
    }
  // find matches and put in final list
  sorted=""
  for (int i = 0; i < order.size(); i++) {
      sorted+=(string_list.find{ it.contains('_'+ order[i] + '.vcf')})+"\n"
    }

  """
  echo "$sorted" > ${sampleID}.list
  """
}
