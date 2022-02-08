def help(){
  println """

Parameter	Default			Description 
--workflow	‘rnaseq’		This is the main parameter that will set the pipeline in motion. 
--pubdir 	"../${workflow}"	The directory that the saved outputs will be stored. 
--organize_by 	‘sample’ 		How to organize the output folder structure. The options are by ‘sample’ or by ‘analysis’ 

cacheDir 	'/projects/omics_share/meta/containers'    This is where the Singularity containers reside 
-w		"/fastscratch/nextflow/${params.workflow}" The directory that all intermediary files and nextflow processes utilize. This needs to be a space with lots of memory. 
--cwd 		System.getProperty("user.dir")             Using Java’s native functions get the current working directory that the script is running. This is required to source the directories of files in the bin folder and other folders. 

--gen_org 	'mouse' (other option 'human')
--extension 	'.fastq.gz' 			 file extension
--read_type 	'PE' (other option 'mouse') 
--sample_folder "/projects/compsci/guglib/tmp_pipeline_defaults/wes_truncated_sequences/pe"  Location that your sequences are to be found
--ref_fa 	Mouse: '/projects/compsci/guglib/tmp_pipeline_defaults/mouse_genome/Mus_musculus.GRCm38.dna.toplevel.fa' 
		Human: '/projects/compsci/refdata/Human/hg38/Index_Files/Bowtie2/Homo_sapiens.GRCh38.dna.toplevel_chr_mod_1_22_MT_X_Y.fa' 
		Reference fasta to be used by various progams

--ref_fa_indices 	'/projects/compsci/refdata/Mouse/mm10/Index_Files/BWA/mm10' 
--filter_trim 		"${params.cwd}/bin/shared/filter_trim.py" 
--min_pct_hq_reads 	0.0 
--read_group_pyfile 	“${params.cwd}/bin/shared/read_group_from_fastq.py" 
--stats_agg 		"${params.cwd}/bin/wes/aggregate_stats.py" 
--mismatch_penalty 	8 
--seed_length 		25 
--rsem_ref_prefix 	Mouse: 'mus_musculus' 
			Human: 'Homo_sapiens.GRCh38.dna.toplevel_chr_mod_1_22_MT_X_Y' 

--rsem_ref_files 	Mouse: '/projects/compsci/guglib/tmp_pipeline_defaults/mouse_rsem_ref' 
			Human: '/projects/compsci/refdata/Human/hg38/Index_Files/Bowtie2' 

--rsem_aligner 		'bowtie2' 
--picard_dict 		Mouse: '/projects/compsci/guglib/tmp_pipeline_defaults/mouse_picard_dict/Mus_musculus.GRCm38.dna.toplevel.dict' 
			Human:'/projects/compsci/refdata/Human/hg38/Index_Files/Bowtie2/Homo_sapiens.GRCh38.dna.toplevel_chr_mod_1_22_MT_X_Y.dict' 
--summary_mets_PE 	"${params.cwd}/bin/rnaseq/summary_QC_metrics_without_xenome.pl" 
--summary_mets_SE 	"${params.cwd}/bin/rnaseq/summary_QC_metrics_without_xenome_SE.pl" 
--probes 		'/projects/compsci/refdata/Human/agilent/hg38_agilent_SureSelect_V4_pChrM_probes_genename.bed' 
--ctp_genes 		'/projects/compsci/refdata/Human/agilent/359genes_b38_noheader_withNames.bed 
--gatk_form 		"${params.cwd}/bin/rnaseq/gatk_formatter.sh" 
--cov_calc 		"${params.cwd}/bin/rnaseq/coveragecalculator.py" 

"""
}
