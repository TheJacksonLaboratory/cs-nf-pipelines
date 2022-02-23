def help(){
  println """

Parameter	Default                 Description
--workflow	‘rnaseq’                This is the main parameter that will set the pipeline in motion.
--pubdir        "../${workflow}"        The directory that the saved outputs will be stored.
--organize_by   ‘sample’                How to organize the output folder structure. The options are by ‘sample’ or by ‘analysis’

cacheDir        '/projects/omics_share/meta/containers'    This is where the Singularity containers reside
-w              "/fastscratch/nextflow/${params.workflow}" The directory that all intermediary files and nextflow processes utilize. 
--cwd           System.getProperty("user.dir")             Using Java’s native functions get the current working directory.

--gen_org	'mouse' (other option 'human')
--extension     '.fastq.gz'                      file extension
--read_type     'PE' (other option 'mouse')
--sample_folder "/projects/compsci/guglib/tmp_pipeline_defaults/wes_truncated_sequences/pe"
--ref_fa        Mouse: '/projects/compsci/guglib/tmp_pipeline_defaults/mouse_genome/Mus_musculus.GRCm38.dna.toplevel.fa'
                Human: '/projects/compsci/refdata/Human/hg38/Index_Files/Bowtie2/Homo_sapiens.GRCh38.dna.toplevel_chr_mod_1_22_MT_X_Y.fa'
                Reference fasta to be used by various progams

--ref_fa_indices        '/projects/compsci/refdata/Mouse/mm10/Index_Files/BWA/mm10'
--filter_trim           "${params.cwd}/bin/shared/filter_trim.py"
--min_pct_hq_reads	0.0
--read_group_pyfile     “${params.cwd}/bin/shared/read_group_from_fastq.py"

--dbSNP 		Mouse: '/projects/compsci/refdata/Mouse/mm10/Annotation_Files/dbSNP/dbSNP.mm10.vcf.gz' 
			Human: '/projects/compsci/refdata/Human/hg38/Annotation_Files/dbsnp/dbsnp_154_hg38.vcf.gz' 
			This points at a database of single nucleotide polymorphisms to be used by multiple GATK programs. 

--target_gatk 		Mouse: '/projects/compsci/refdata/Mouse/mm10/mm10Exome_v4_12-19.1.mm10.baits_merged.bed' 
			Human: '/projects/compsci/refdata/Human/agilent/hg38_liftover_agilent_SureSelect_V4_pChrM_probes.bed' 
			A bed file that is used for targeting specific genes when working with GATK 

--target_picard 	Mouse: '/projects/compsci/refdata/Mouse/mm10/mm10Exome_v4_12-19.1.mm10.targets_merged_picard_new.bed' 
			Human: '/projects/compsci/refdata/Human/agilent/hg38_agilent_SureSelect_V4_pChrM_probes_picard_updated.bed' 
			a bed file used when targeting genes using Picard tools 

--bait_picard 		Mouse: '/projects/compsci/refdata/Mouse/mm10/mm10Exome_v4_12-19.1.mm10.baits_merged_picard_new.bed' 
			Human: '/projects/compsci/refdata/Human/agilent/hg38_agilent_SureSelect_V4_pChrM_probes_picard_updated.bed' 
			A bed file used to bait jeans while using Picard 

--stats_agg 		"${params.cwd}/bin/wes/aggregate_stats.py" 
--mismatch_penalty 	8 
--call_val 		50 	The minimum phred-scaled confidence threshold at which variants should be called. 
--ploidy_val 		Mouse: Null 
			Human: "-ploidy 4" 
			Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy 
--gen_ver 		"hg38" 
--snpEff_config 	"/projects/compsci/refdata/Mouse/mm10/snpEff_files/snpEff.config" 	The configuration file used while running snpEff 
			NOTE: snpEff.config needs to be in directory one level above data directory containing snpEff database files 
--gold_std_indels 	'/projects/compsci/refdata/Human/hg38/Annotation_Files/indels/Mills_and_1000G_gold_standard.indels.hg38.vcf’ 
--phase1_1000G 		'/projects/compsci/refdata/Human/hg38/Annotation_Files/indels/1000G_phase1.snps.high_confidence.hg38.vcf' 
--dbNSFP 		'/projects/compsci/refdata/Human/hg38/Annotation_Files/dbNSFP/hg38_dbNSFP3.2a.txt.gz' 
--cosmic 		'/projects/compsci/refdata/Human/hg38/Annotation_Files/Cosmic/COSMICv92_Coding_NonCoding.vcf' 
--cosmic_annot 		“${params.cwd} /bin/Cosmic_Annotation_hg38.pl” 
--hgvs_data 		'/projects/compsci/refdata/Human/hg38/snpEff_files/data/hg38' 
"""
}

