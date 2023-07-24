def help(){
  println '''
Parameter | Type | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--organize_by | sample | How to organize the output folder structure. Options: sample or analysis.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /fastscratch or other directory with ample storage.

--sample_folder | /<PATH> | The path to the folder that contains all the samples to be run by the pipeline. The files in this path can also be symbolic links. 
--extension | <string> | Default: '.fastq.gz' The expected extension for the input read files.
--pattern | <string> | Default: '*_R{1,2}*'. The expected R1 / R2 matching pattern. The default value will match reads with names like this READ_NAME_R1_MoreText.fastq.gz or READ_NAME_R1.fastq.gz
--read_type | <string> | Options: PE and SE. Default: PE. Type of reads: paired end (PE) or single end (SE).
--csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.
--download_data | null | Requires `--csv_input`. When specified, read data in the CSV manifest will be downloaded from provided URLs. 
--concat_lanes | <boolean> | Options: false and true. Default: false. If this boolean is specific, FASTQ files will be concatenated by sample. This option is used in cases where samples are divided across individual sequencing lanes.

NOTE: When `--concat_lanes` is used. Unique Sample IDs must be parsed from FASTQ names. The following commands split FASTQ names on a delimiter and keep 'n' positions of the split array. 

--concat_sampleID_delim | <string> | Default: '_'. The delimited to split FASTQ file names. 
--concat_sampleID_positions | <numeric> | Default: 1. The number of elements to keep after splitting on the chosen delimiter in the sample name. 

Examples: 
    Given the input file name "SAMPLE_NAME_1_OTHER_STUFF-WeDont_WANT.txt" if this `concat_sampleID_delim` = '_' and `concat_sampleID_positions` = "3" the sample ID would be assigned as `SAMPLE_NAME_1`
    Given the input file name "SAMPLE_NAME_1_OTHER_STUFF-WeDont_WANT.txt" if this `concat_sampleID_delim` = '-' and `concat_sampleID_positions` = "1" the sample ID would be assigned as `SAMPLE_NAME_1_OTHER_STUFF`

--bowtie_index | /<PATH> | Path to the bowtie index. Include the bowtie prefix in this path (e.g., `/path/to/bowtie.transcripts` where bowtie.transcripts.* are the full set of index files in the directory.  
--transcripts_info | /<PATH> | A file containing all transcript IDs. NOTE: These IDs must not contain haplotype IDs. This file must also have a 'length' column. Note that 'length' is not used in this context. ONLY IDs are used from this file. Can be obtained from `prepare_emase` workflow (emase.fullTranscripts.info)
--gbrs_strain_list | <comma,delim,list> | A list of haplotype names corresponding to genomes used in hybrid genome construction (e.g., 'A,B,C,D,E,F,G,H'). 
--gene2transcript_csv | /<PATH> | A file containing all gene to transcript ID translations. NOTE: These IDs must not contain haplotype IDs. Can be obtained from `prepare_emase` workflow (emase.gene2transcripts.tsv)
--full_transcript_info | /<PATH> | A file containing all transcript IDs with transcript lengths. NOTE: These IDs must contain haplotype IDs. Can be obtained from `prepare_emase` workflow (emase.pooled.fullTranscripts.info)
--emase_model | <numeric> | Options: 1, 2, 3, 4. 
  1: reads are apportioned among genes first,
      then between alleles, and then among isoforms.
  2: reads are apportioned among genes first,
      then among isoforms, and then between alleles.
  3: reads are apportioned among genes first,
      then among each isoform-allele combination
      which are treated equally.
  4: assumes no hierarchy and multi-reads are
      apportioned equally among genes, isoforms, and
      alleles

--keep_intermediate | <boolean> | Default: false. Keep intermediate files, not otherwise saved. 
'''
}
