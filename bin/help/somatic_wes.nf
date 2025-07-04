def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--sample_folder | /<PATH> | The path to the folder that contains all the samples to be run by the pipeline. The files in this path can also be symbolic links. 
--extension | .fastq.gz | The expected extension for the input read files.
--pattern | '*_R{1,2}*' | The expected R1 / R2 matching pattern. The default value will match reads with names like this READ_NAME_R1_MoreText.fastq.gz or READ_NAME_R1.fastq.gz
--read_type | PE | Options: PE and SE. Default: PE. Type of reads: paired end (PE) or single end (SE).
--concat_lanes | false | Options: false and true. Default: false. If this boolean is specified, FASTQ files will be concatenated by sample. This option is used in cases where samples are divided across individual sequencing lanes.
--csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.
--download_data | null | Requires `--csv_input`. When specified, read data in the CSV manifest will be downloaded from provided URLs. 

--gen_org | human | Options: human only.

--quality_phred | 15 | The quality value that is required for a base to pass. Default: 15 which is a phred quality score of >=Q15.
--unqualified_perc | 40 | Percent of bases that are allowed to be unqualified (0~100). Default: 40 which is 40%.
--detect_adapter_for_pe | false | If true, adapter auto-detection is used for paired end data. By default, paired-end data adapter sequence auto-detection is disabled as the adapters can be trimmed by overlap analysis. However, --detect_adapter_for_pe will enable it. Fastp will run a little slower if you specify the sequence adapters or enable adapter auto-detection, but usually result in a slightly cleaner output, since the overlap analysis may fail due to sequencing errors or adapter dimers.

--pdx | false | Options: false, true. If specified, 'Xengsort' is run on reads to deconvolute human and mouse reads. Human only reads are used in analysis. 
--xengsort_host_fasta | '/projects/compsci/omics_share/mouse/GRCm39/genome/sequence/imputed/rel_2112_v8/NOD_ShiLtJ.39.fa' | Xengsort host fasta file. Used by Xengsort Index when `--pdx` is run, and xengsort_idx_path is `null` or false.  
--xengsort_idx_path = | '/projects/compsci/omics_share/human/GRCh38/supporting_files/xengsort' | Xengsort index for deconvolution of human and mouse reads. Used when `--pdx` is run. If `null`, Xengsort Index is run using ref_fa and host_fa.  
--xengsort_idx_name = | 'hg38_GRCm39-NOD_ShiLtJ' | Xengsort index name associated with files located in `xengsort_idx_path` or name given to outputs produced by Xengsort Index

--genotype_targets | '/projects/compsci/omics_share/human/GRCh38/supporting_files/ancestry_panel/snp_panel_v2_targets_annotations.snpwt.bed.gz' | Target SNP bed file for the ancestry panel. Can contain annotation information. 
--snpID_list | '/projects/compsci/omics_share/human/GRCh38/supporting_files/ancestry_panel/snp_panel_v2.list' | Target SNPs in list used in BCFtools filtering step
--snp_annotations | '/projects/compsci/omics_share/human/GRCh38/supporting_files/ancestry_panel/snp_panel_v2_targets_annotations.snpwt.bed.gz' | Target SNP bed file with annotations for the ancestry panel.
--snpweights_panel | '/projects/compsci/omics_share/human/GRCh38/supporting_files/ancestry_panel/ancestry_panel_v2.snpwt' | SNP weights panel in the appropriate format

--ref_fa | '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta' | The reference fasta to be used throughout the process for alignment as well as any downstream analysis. JAX users should not change this parameter.
--ref_fa_indices | '/projects/omics_share/human/GRCh38/genome/indices/gatk/bwa/Homo_sapiens_assembly38.fasta' | Pre-compiled BWA index files. JAX users should not change this parameter.
--target_gatk | '/projects/omics_share/human/GRCh38/supporting_files/capture_kit_files/agilent/v7/S31285117_MergedProbes_no_gene_names.bed' | A bed file with WES target intervals as defined in the capture array used in the data. NOTE: This file MUST reflect the capture array used to generate your data.
--target_picard | '/projects/omics_share/human/GRCh38/supporting_files/capture_kit_files/agilent/v7/S31285117_MergedProbes_no_gene_names.picard.interval_list' | A GATK interval file covering WES target intervals. Used in calculating coverage metrics. NOTE: This file MUST reflect the capture array used to generate your data.
--bait_picard | '/projects/omics_share/human/GRCh38/supporting_files/capture_kit_files/agilent/v7/S31285117_MergedProbes_no_gene_names.picard.interval_list' | A GATK interval file covering WES target intervals. Used in calculating coverage metrics. This file can be the same as the interval file,  NOTE: This file MUST reflect the capture array used to generate your data.

--mismatch_penalty | 8 | The BWA penalty for a mismatch.

--gnomad_ref | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/snps_indels/af-only-gnomad.hg38.vcf.gz' | GnomAD germline reference from GATK resource pack. 
--pon_ref | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/snps_indels/1000g_pon.hg38.vcf.gz' | 1000 genome germline panel of normals from GATK resource pack.
--genotype_pon | true | Call sites in the PoN even though they will ultimately be filtered.
--genotype_germline | true | Call all apparent germline site even though they will ultimately be filtered.
--contam_ref | '/projects/compsci/omics_share/human/GRCh38/genome/annotation/snps_indels/small_exac_common_3.hg38.vcf.gz' | File used in GetPileupSummaries and CalculateContaminationcommon. A germline variant sites VCF, e.g. derived from the gnomAD resource, with population allele frequencies (AF) in the INFO field is used from GATK resource bundle.

--dbSNP | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz' | The dbSNP database contains known single nucleotide polymorphisms, and is used in the annotation of known variants. Points to human dbSNP when --gen_org human. JAX users should not change this parameter.

--gen_ver | 'hg38' | snpEff genome version. Sets to 'hg38' when --gen_org human JAX users should not change this parameter.

--snpEff_config | Human: '/projects/omics_share/human/GRCh38/genome/indices/snpEff_5_1/snpEff.config' | The configuration file used while running snpEff, points to human snpEff file when --gen_org human. JAX users should not change this parameter.

--gold_std_indels | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz’ | Human Only - Used in GATK BaseRecalibrator. JAX users should not change this parameter.
--phase1_1000G | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/1000G_phase1.snps.high_confidence.hg38.vcf.gz' | Human Only - Used in GATK BaseRecalibrator. JAX users should not change this parameter.
--dbNSFP | '/projects/omics_share/human/GRCh38/genome/annotation/function/dbNSFP4.2a.gatk_formatted.txt.gz' | Human Only - Used in variant annotation.
--cosmic | '/projects/omics_share/human/GRCh38/genome/annotation/function/COSMICv95_Coding_Noncoding.gatk_formatted.vcf' | Human Only - Used in variant annotation.
'''
}

