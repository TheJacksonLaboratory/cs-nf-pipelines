def help(){
  println '''
Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /fastscratch or other directory with ample storage.
--snp_vcf | /<PATH> | VCF containing only SNPs for patching into primary reference fasta. 
--indel_vcf | /<PATH> | VCF containing only InDELs for transforming into primary reference fasta. 
--primary_reference_fasta | /<PATH> | The primary reference fasta file, where patched SNPs and transformed InDELs are applied.
--primary_reference_gtf | /<PATH> | The primary reference gtf file, used to patch and transform gene/transcripts/exons. 
--gtf_biotype_include | protein_coding,lncRNA,IG_C_gene,IG_D_gene,IG_J_gene,IG_LV_gene,IG_V_gene,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene | A comma delimited list of terms to include from the full GTF. All other biotype terms will be excluded. 
--append_chromosomes | true | Add back any full chromosomes that are dropped due to lack of variants in the SNP or INDEL file. See the Wiki for additional information.
--strain | <comma,delim,string> | A comma delimited string of strains/haplotypes. (e.g., 'A_J,CAST_EiJ,...')
--genome_version | <string> | A genome ID string (e.g., 39)
--diploid | <boolean> | Create diploid VCI file
--keep_fails | <boolean> | Default: false. Keep track of VCF lines that could not be converted to VCI file
--pass_only | <boolean> | Default: false. Use only VCF lines that have a PASS for the filter value
--quality_filter | <string> | Default: NULL. Filter on quality, (e.g., 'FI=PASS')
--region | <seqid:start-end> | Default: NULL. A region used in extraction. If using this option, the bed option can not be used.
--bed | /<BED_PATH> | Default: NULL. A BED file with regions for extraction. This option cannot be used with region.

--keep_intermediate | <boolean> | Default: false. Keep intermediate files, not otherwise saved. 
'''
}
