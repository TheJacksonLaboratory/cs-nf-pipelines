process MAKE_VCF_LIST {

 
 input:
 tuple val(sampleID), val(chromes)

 output:
 tuple val(sampleID), file("*.list"), emit: list

 script:
 log.info "----- Make VCF List from Chromosomes: ${sampleID} -------"
 def my_list=chromes.join("\n")
 """
 echo "$my_list" > my.list
 """
}
