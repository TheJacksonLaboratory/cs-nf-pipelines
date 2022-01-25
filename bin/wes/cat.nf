process CAT_SNP_INDEL{
  """
  cat !{snps_filt} > !{sampleID}_full_anno_snp.vcf
  cat !{indels_filt} > !{sampleID}_full_anno_indel.vcf
  """
}
process CAT_INDEL_HUMAN{
 """
 cat !{sampleID}_full_anno_indel.vcf | /snpEff_v4_3/snpEff/scripts/vcfEffOnePerLine.pl > !{sampleID}_full_anno_indel_onePerLine.vcf
"""
}
