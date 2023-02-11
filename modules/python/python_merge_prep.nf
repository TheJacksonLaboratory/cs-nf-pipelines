process MERGE_PREP {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), path(vcf), val(meta), val(normal_name), val(tumor_name), val(tool)
  
  output:
  tuple val(sampleID), path("*_mergePrep.vcf"), val(meta), val(normal_name), val(tumor_name), val(tool), emit: merge_prep_vcf

  script:
  String support_call = tool == 'manta' || tool == 'lancet_support' ? '--support' : ''
  String tool_name = tool == 'lancet_support' ? 'lancet' : tool

  """
  python \
  ${projectDir}/bin/sv/reorder_vcf.py \
  ${vcf} \
  ${vcf.baseName}_ordered.vcf \
  ${normal_name} ${tumor_name}
  
  python \
  ${projectDir}/bin/sv/merge_prep.py \
  --vcf ${vcf.baseName}_ordered.vcf \
  --out ${vcf.baseName}_mergePrep.vcf \
  --tool ${tool_name} \
  ${support_call}
  """
}

/* NOTE: PLEASE READ!!! 
  `reorder_vcf.py` requires the header and input 'tumor/normal' names in the 3rd and 4th arg to match. 
      If you pass names NOT present in the header, it will simply emit the file AS IS.
      The script DOES NOT inform the user of if a change has been made in the sample order. 
      NOTE ALSO: if the header already contains the strings 'TUMOR' and 'NORMAL,
                 'TUMOR and NORMAL are RENAMED to string provided in 3rd and 4th args. 
*/
