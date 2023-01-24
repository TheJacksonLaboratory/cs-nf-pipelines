def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /fastscratch or other directory with ample storage.

--genome_file_list | /<PATH> OR /<PATH>,/<PATH/,... | A comma seperated list of FASTA genome file(s) for use hybrid genome construction (e.g., genome1.fa OR genome1.fa,genome2.fa,...). NOTE: FASTA AND GTF MUST BE IN THE SAME ORDER.  
--gtf_file_list | /<PATH> OR /<PATH>,/<PATH/,... |  A comma seperated list of GTF files corresponding to the genomes for use hybrid transcriptome construction (e.g., genome1.gtf OR genome1.gtf,genome2.gtf,...). NOTE: GTF AND FASTA MUST BE IN THE SAME ORDER. 
--haplotype_list | <comma,delim,string> | A list of haplotype names corresponding to genomes used in hybrid genome contrucution (e.g., 'A,B,C,D,E,F,G,H'). These names are appended to transcript IDs (e.g., ENMST00000042_A). NOTE: HAPLOTYPE LIST MUST BE IN THE SAME ORDER AS FASTA AND GTF FILES. 

--keep_intermediate               ${params.keep_intermediate}

'''
}
