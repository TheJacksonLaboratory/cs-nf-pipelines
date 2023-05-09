def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--organize_by | sample | How to organize the output folder structure. Options: sample or analysis.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /fastscratch or other directory with ample storage.

--sample_folder | /<PATH> | The path to the folder that contains all the samples to be run by the pipeline. The files in this path can also be symbolic links. 
--extension | .fastq.gz | The expected extension for the input read files.
--pattern | '*_R{1,2}*' | The expected R1 / R2 matching pattern. The default value will match reads with names like this READ_NAME_R1_MoreText.fastq.gz or READ_NAME_R1.fastq.gz
--read_type | PE | Options: PE and SE. Default: PE. Type of reads: paired end (PE) or single end (SE).
--concat_lanes | false | Options: false and true. Default: false. If this boolean is specified, FASTQ files will be concatenated by sample. This option is used in cases where samples are divided across individual sequencing lanes.

--gen_org | mouse | Options: mouse and human.

--effective_genome_size | The length of the “mappable” genome. | Mouse only - Please see : 'https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html'.
                        
--chain | null | The default value for Mouse Reference Strain -  g2gtools chain file to adjust coordinates to reference                                                              

--bowtie2Index | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/ensembl/v102/bowtie2/Mus_musculus.GRCm38.dna.primary_assembly.fa' | The default value for mm10. 
              | Human: '/projects/omics_share/human/GRCh38/genome/indices/gatk/bowtie2/hg38_noalt'.
              | Pre-compiled BOWTIE2 index files, points to human reference when --gen_org human.

--bowtieMaxInsert | 1000 |  The maximum fragment length for valid paired-end alignments.

--cutadaptMinLength | 20 | The minimum length to discard processed reads.
--cutadaptQualCutoff | 20 | The quality cutoff used to trim low-quality ends from reads.
--cutadaptAdapterR1 | Specification of a 3’ adapter or a linked adapter.
--cutadaptAdapterR2 | Specification of a 5’ adapter or a linked adapter.
--tmpdir  | /<PATH> | Temporary directory to store temp files.  
'''
}

