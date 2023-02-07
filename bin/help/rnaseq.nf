def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--organize_by | sample | How to organize the output folder structure. Options: sample or analysis
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /fastscratch or other directory with ample storage.

--sample_folder | /<PATH> | The path to the folder that contains all the samples to be run by the pipeline. The files in this path can also be symbolic links. 
--extension | .fastq.gz | The expected extension for the input read files.
--pattern | '*_R{1,2}*' | The expected R1 / R2 matching pattern. The default value will match reads with names like this READ_NAME_R1_MoreText.fastq.gz or READ_NAME_R1.fastq.gz
--read_type | PE | Options: PE and SE. Default: PE. Type of reads: paired end (PE) or single end (SE).
--concat_lanes | false | Options: false and true. Default: false. If this boolean is specific, FASTQ files will be concatenated by sample. This option is used in cases where samples are divided across individual sequencing lanes.

--gen_org | mouse | Options: mouse and human.

--pdx | false | Options: true or false. If 'true' Xenome is run to remove mouse reads from samples. 
--xenome_prefix | '/projects/compsci/omics_share/human/GRCh38/supporting_files/xenome/trans_human_GRCh38_84_NOD_based_on_mm10_k25' | Pre-compiled Xenome classification index files. Used if PDX analysis is specified. 

--read_prep | 'reverse_stranded' | Options: 'reverse_stranded', 'forward_stranded' or 'non_stranded'. This determines how RNA quantification is done, and statistics are calculated. It is based on the library strandedness. 

--min_pct_hq_reads| '0.0' | The minimum percent of high-quality reads passing when trimming the fastq files.
--rsem_ref_files | /projects/omics_share/mouse/GRCm38/transcriptome/indices/ensembl/v102/bowtie2 | Pre-compiled index files. Refers to human indices when --gen_org human. JAX users should not change this, unless using STAR indices.
--rsem_ref_prefix | 'Mus_musculus.GRCm38.dna.toplevel' | Prefix for index files. JAX users should not change this, unless using STAR indices. Refers to human indices when --gen_org human.
--seed_length | 25 | 'Seed length used by the read aligner. Providing the correct value is important for RSEM. If RSEM runs Bowtie, it uses this value for Bowtie's seed length parameter.'
--rsem_aligner | 'bowtie2' | Options: bowtie2 or star. The aligner algorithm used by RSEM. Note, if using STAR, point rsem_ref_files to STAR based indices.

--picard_dict | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.dict' 
              | Human: '/projects/omics_share/human/GRCh38/genome/sequence/ensembl/v104/Homo_sapiens.GRCh38.dna.toplevel.dict'
              | The coverage metric calculation step requires this file. Refers to human assembly when --gen_org human. JAX users should not change this parameter.

--ref_flat | Mouse: '/projects/omics_share/mouse/GRCm38/transcriptome/annotation/ensembl/v102/Mus_musculus.GRCm38.102.chr_patch_hapl_scaff.refFlat.txt' 
           | Human: '/projects/omics_share/human/GRCh38/transcriptome/annotation/ensembl/v104/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.refFlat.txt'
           | The coverage metric calculation step requires this file. Refers to human assembly when --gen_org human. JAX users should not change this parameter.

--ribo_intervals | Mouse: '/projects/omics_share/mouse/GRCm38/transcriptome/annotation/ensembl/v102/Mus_musculus.GRCm38.102.chr_patch_hapl_scaff.rRNA.interval_list' 
                 | Human: '/projects/omics_share/human/GRCh38/transcriptome/annotation/ensembl/v104/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.rRNA.interval_list'
                 | The coverage metric calculation step requires this file. Refers to human assembly when --gen_org human. JAX users should not change this parameter.

There are two additional parameters that are human specific. They are: 

Parameter| Default| Description

--ref_fa | '/projects/omics_share/human/GRCh38/genome/sequence/ensembl/v104/Homo_sapiens.GRCh38.dna.toplevel.fa'| Reference fasta to be used in alignment calculation as well as any downstream analysis. JAX users should not change this parameter.
--ref_fai | '/projects/omics_share/human/GRCh38/genome/sequence/ensembl/v104/Homo_sapiens.GRCh38.dna.toplevel.fa.fai' | Reference fasta index file.  JAX users should not change this parameter.
'''
}
