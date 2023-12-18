def help(){
  println '''
--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--num_generations | 100 | The number of generations to calculate transition probabilities out to.  
--ensembl_build | 105 | The ensembl build number used to extract gene names and locations from the R package `biomaRt`.
--emase_gene2transcript | /<PATH> | A file containing all gene to transcript ID translations. NOTE: These IDs must not contain haplotype IDs. Can be obtained from `prepare_emase` workflow (emase.gene2transcripts.tsv)
'''
}
