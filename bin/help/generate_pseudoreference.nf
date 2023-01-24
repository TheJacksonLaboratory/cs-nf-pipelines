def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /fastscratch or other directory with ample storage.
--snp_vcf | /<PATH> | VCF containing only SNPs for patching into primary reference fasta. 
--indel_vcf | /<PATH> | VCF containing only InDELs for transforming into primary reference fasta. 
--primary_reference_fasta | /<PATH> | The primary reference fasta file, where patched SNPs and transformed InDELs are applied.
--primary_reference_gtf | /<PATH> | The primary reference gtf file, used to patch and transform gene/transcripts/exons. 
--strain | <comma,delim,string> | A comma delimilited string of strains/haplotypes. (e.g., 'A_J,CAST_EiJ,...')
--genome_version | <string> | A genome ID string (e.g., 39)
--diploid | boolean | Create diploid VCI file
--keep_fails | boolean | Keep track of VCF lines that could not be converted to VCI file
--pass_only | boolean | Use only VCF lines that have a PASS for the filter value
--quality_filter | <string> | Filter on quality, (e.g., 'FI=PASS')
--region | <seqid:start-end> | A region used in extraction. If using this option, the bed option can not be used.
--bed | /<BED_PATH> | A BED file with regions for extraction. This option cannot be used with region.

'''
}
