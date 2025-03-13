def help(){
  println '''
Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--ref_fa | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.fa' or GRCm39: '/projects/omics_share/mouse/GRCm39/genome/sequence/ensembl/v105/Mus_musculus.GRCm39.dna.primary_assembly.fa
         | Human: '/projects/omics_share/human/GRCh38/genome/sequence/ensembl/v104/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
         | The reference fasta to be used in index generation.

--ref_gtf | Mouse: '/projects/omics_share/mouse/GRCm38/transcriptome/annotation/ensembl/v102/Mus_musculus.GRCm38.102.gtf' or GRCm39: '/projects/omics_share/mouse/GRCm39/transcriptome/annotation/ensembl/v105/Mus_musculus.GRCm39.105.gtf
          | Human: '/projects/omics_share/human/GRCh38/transcriptome/annotation/ensembl/v104/Homo_sapiens.GRCh38.104.gtf'
          | The reference fasta to be used in index generation.

--custom_gene_fasta | null | The path to a fasta file with additonal transcript sequences to add to the index. Will be annotated based on the name provided in the sequnece name field. For example: ">New_Gene_42", where New_Gene_42 will be the name of the gene, transcript, and exon. 

'''
}
