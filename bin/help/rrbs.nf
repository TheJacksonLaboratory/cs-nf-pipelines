def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.                                                                                                                                               |
--organize_by | sample | How to organize the output folder structure. Options: sample or analysis                                                                                                                       |
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.                                                                                         |
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /fastscratch or other directory with ample storage. |

--sample_folder | /<PATH> | The path to the folder that contains all the samples to be run by the pipeline. The files in this path can also be symbolic links. 
--extension | .fastq.gz | The expected extension for the input read files.
--pattern | '*_R{1,2}*' | The expected R1 / R2 matching pattern. The default value will match reads with names like this READ_NAME_R1_MoreText.fastq.gz or READ_NAME_R1.fastq.gz
--read_type | PE | Options: PE and SE. Default: PE. Type of reads: paired end (PE) or single end (SE).
--concat_lanes | false | Options: false and true. Default: false. If this boolean is specified, FASTQ files will be concatenated by sample. This option is used in cases where samples are divided across individual sequencing lanes.
--csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.
--download_data | null | Requires `--csv_input`. When specified, read data in the CSV manifest will be downloaded from provided URLs. 

--gen_org | mouse | Options: mouse and human.
--genome_build | 'GRCm38' | Mouse specific. Options: GRCm38 or GRCm39. If gen_org == human, build defaults to GRCm38.

--non_directional | true | Options: true and false. Selecting this option for non-directional RRBS libraries will screen quality-trimmed sequences for CAA or CGA at the start of the read and, if found, removes the first two base pairs. 

--trimLength | 30 | Discard reads that became shorter than length 'INT' because of either quality or adapter trimming. A value of 0 effectively disables this behaviour.
--qualThreshold | 30 | Trim low-quality ends from reads in addition to adapter removal. For RRBS samples, quality trimming will be performed first, and adapter trimming is carried in a second round. Other files are quality and adapter trimmed in a single pass. The algorithm is the same as the one used by BWA (Subtract INT from all qualities; compute partial sums from all indices to the end of the sequence; cut sequence at the index at which the sum is minimal).
--adapOverlap | 1 | Stringency for overlap with adapter sequence required to trim a sequence. Defaults to a very stringent setting of 1, i.e. even a single base pair of overlapping sequence will be trimmed of the 3' end of any read.
--adaptorSeq | 'AGATCGGAAGAGC' | Adapter sequence to be trimmed. This sequence is the standard Illumina adapter sequence. 

--seedLength | 20 | Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive    
--seedMismatch | 0 | Sets the number of mismatches to be allowed in a seed alignment during multiseed alignment. Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower) but increases sensitivity.
--MinInsert | 0 | The minimum insert size for valid paired-end alignments. E.g. if -I 60 is specified and a paired-end alignment consists of two 20-bp alignments in the appropriate orientation with a 20-bp gap between them, that alignment is considered valid (as long as -X is also satisfied). A 19-bp gap would not be valid in that case.
--MaxInsert | 1000 | The maximum insert size for valid paired-end alignments. E.g. if -X 100 is specified and a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied). A 61-bp gap would not be valid in that case

--ref_fa_index | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/ensembl/v102/bismark/bowtie2'
               | Human: '/projects/omics_share/human/GRCh38/genome/indices/gatk/bismark/bowtie2'
               | Pre-compiled Bismark Bowtie2 index files. Points to human reference when --gen_org human. JAX users should change this parameter only if changing aligner used.
--aligner | 'bowtie2' | Options: bowtie2 and hisat2. Bismark alignment tool. If hisat2 is specified, change `ref_fa_index` to the appropriate index files. 

--skip_deduplication | true | Bismark Deduplication is used if `true`

--cytosine_report | false | After the conversion to bedGraph has completed, the option --cytosine_report produces a genome-wide methylation report for all cytosines in the genome
--comprehensive | true | Specifying this option will merge all four possible strand-specific methylation info into context-dependent output files.
'''

}