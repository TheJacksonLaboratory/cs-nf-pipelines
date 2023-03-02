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
--concat_lanes | false | Options: false and true. Default: false. If this boolean is specific, FASTQ files will be concatenated by sample. This option is used in cases where samples are divided across individual sequencing lanes.
--csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.
--download_data | null | Requires `--csv_input`. When specified, read data in the CSV manifest will be downloaded from provided URLs. 

--gen_org | mouse | Options: mouse and human.

--xenome_prefix | /projects/compsci/omics_share/human/GRCh38/supporting_files/xenome/trans_human_GRCh38_84_NOD_based_on_mm10_k25| Xenome index for deconvolution of human and mouse reads. Used when `--pdx` is run. 
--read_length | 150 | Options: 75, 100, 150. Changed relative to sample read length.
--star_index | /projects/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/star/star-2.7.4a-150bp | STAR index used by several tools. Change the index relative to sample read length. Read length options: 75, 100, 150. 
--star_fusion_star_index | /projects/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/starfusion/star-150 | STAR-fusion index. Change the index relative to sample read length. Read length options: 75, 100, 150. 

--gtf | /projects/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/ensembl/Homo_sapiens.GRCh38.102.gtf | GTF file used by several callers. STAR refrences were built from this file, and it should not be changed unless other indicies are also updated. 
--fasta | /projects/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/ensembl/Homo_sapiens.GRCh38.102.all.fa | Genomic FASTA file used by fusion callers. STAR refrences were built from this file, and it should not be changed unless other indicies are also updated. 

--arriba_star_args | <see_config_file> | Arriba recommended argument string for STAR alignment. See the rna_fusion.config file for specific arguments used. 
--arriba_blacklist | /projects/compsci/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/arriba/blacklist_hg38_GRCh38_v2.4.0.tsv.gz | Arriba provided blacklist of difficult regions for fusion calling: https://arriba.readthedocs.io/en/latest/input-files/
--arriba_known_fusions | /projects/compsci/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/arriba/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz | Arriba provided list of known fusions: https://arriba.readthedocs.io/en/latest/input-files/
--arriba_protein_domains | /projects/compsci/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3 | Arriba provided known protein domains: https://arriba.readthedocs.io/en/latest/input-files/

--fusioncatcher_ref | /projects/compsci/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/fusioncatcher/human_v102 | Fusion catcher provided reference files: http://sourceforge.net/projects/fusioncatcher/files/data/
--fusioncatcher_limitSjdbInsertNsj | 2000000 | STAR option used by fusioncatcher: maximum number of junction to be inserted to the genome on the fly at the mapping stage 

--jaffa_ref_dir | /projects/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/jaffa/ | Jaffa provided reference files: https://github.com/Oshlack/JAFFA/wiki/Download

--kallisto_index | /projects/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/pizzly/Homo_sapiens.GRCh38.102.cdna.all.kallisto-0.48.0.index | Kallisto alignment index. 
--transcript_fasta | /projects/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/ensembl/Homo_sapiens.GRCh38.102.cdna.all.fa.gz | Transcriptome FASTA file used by Pizzly. 

--squid_star_args | <see_config_file> | Squid recommended argument string for STAR alignment. See the rna_fusion.config file for specific arguments used. 

--star_fusion_ref | /projects/omics_share/human/GRCh38/transcriptome/indices/rna_fusion/starfusion/ctat_genome_lib_build_dir | star-fusion reference file set. Build from the above GTF and FASTA.
--star_fusion_opt | null | Additional star-fusion options can be provided. 

--fusion_report_opt | null | Additional fusion-report options can be provided. 
--databases | /projects/compsci/omics_share/human/GRCh38/supporting_files/rna_fusion_dbs | Fusion-report databases of known fusion events. Used in report generation only. 

--pdx | false | Options: false, true. If specified, 'Xenome' is run on reads to deconvolute human and mouse reads. Human only reads are used in analysis. 

'''
}

