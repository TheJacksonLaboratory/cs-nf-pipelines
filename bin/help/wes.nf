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

--gen_org | mouse | Options: mouse and human.

--ref_fa | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.toplevel.fa' 
         | Human: '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta'
         | The reference fasta to be used throughout the process for alignment as well as any downstream analysis, points to human reference when --gen_org human. JAX users should not change this parameter.

--ref_fa_indices | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/ensembl/v102/bwa/Mus_musculus.GRCm38.dna.toplevel.fa' 
                 | Human: '/projects/omics_share/human/GRCh38/genome/indices/gatk/bwa/Homo_sapiens_assembly38.fasta.64'
                 | Pre-compiled BWA index files, points to human reference when --gen_org human. JAX users should not change this parameter.

--min_pct_hq_reads | 0.0 | The minimum percent of high-quality reads passing when trimming the fastq files.

--target_gatk | Mouse: '/projects/omics_share/mouse/GRCm38/supporting_files/capture_kit_files/agilent/v2/S32371113_mouse_exon_V2.bare.bed' 
              | Human: '/projects/omics_share/human/GRCh38/supporting_files/capture_kit_files/agilent/v7/S31285117_MergedProbes_no_gene_names.bed'
              | A bed file with WES target intervals as defined in the capture array used in the data, points to human bed file when --gen_org human. NOTE: This file MUST reflect the capture array used to generate your data.

--target_picard | Mouse: '/projects/omics_share/mouse/GRCm38/supporting_files/capture_kit_files/agilent/v2/S32371113_mouse_exon_V2.bare.picard.interval_list' 
                | Human: '/projects/omics_share/human/GRCh38/supporting_files/capture_kit_files/agilent/v7/S31285117_MergedProbes_no_gene_names.picard.interval_list'
                | A GATK interval file covering WES target intervals. Used in calculating coverage metrics. NOTE: This file MUST reflect the capture array used to generate your data.

--bait_picard | Mouse: '/projects/omics_share/mouse/GRCm38/supporting_files/capture_kit_files/agilent/v2/S32371113_mouse_exon_V2.bare.picard.interval_list' 
              | Human: '/projects/omics_share/human/GRCh38/supporting_files/capture_kit_files/agilent/v7/S31285117_MergedProbes_no_gene_names.picard.interval_list'
              | A GATK interval file covering WES target intervals. Used in calculating coverage metrics. This file can be the same as the interval file,  NOTE: This file MUST reflect the capture array used to generate your data.

--mismatch_penalty | -B 8 | The BWA penalty for a mismatch.
--call_val | 50 | The minimum phred-scaled confidence threshold at which variants should be called.
--ploidy_val | '-ploidy 2' | Sample ploidy

--dbSNP | Mouse: '/projects/omics_share/mouse/GRCm38/genome/annotation/snps_indels/GCA_000001635.6_current_ids.vcf.gz' 
        | Human: '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz'
        | The dbSNP database contains known single nucleotide polymorphisms, and is used in the annotation of known variants. Points to human dbSNP when --gen_org human. JAX users should not change this parameter.

--gen_ver | Mouse: 'GRCm38.99' 
          | Human: 'hg38'
          | snpEff genome version. Sets to 'hg38' when --gen_org human JAX users should not change this parameter.

--snpEff_config | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/snpEff_5_1/snpEff.config' 
                | Human: '/projects/omics_share/human/GRCh38/genome/indices/snpEff_5_1/snpEff.config'
                | The configuration file used while running snpEff, points to human snpEff file when --gen_org human. JAX users should not change this parameter.


--gold_std_indels | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz’ | Human Only - Used in GATK BaseRecalibrator. JAX users should not change this parameter.
--phase1_1000G | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/1000G_phase1.snps.high_confidence.hg38.vcf.gz' | Human Only - Used in GATK BaseRecalibrator. JAX users should not change this parameter.
--dbNSFP | '/projects/omics_share/human/GRCh38/genome/annotation/function/dbNSFP4.2a.gatk_formatted.txt.gz' | Human Only - Used in variant annotation.
--cosmic | '/projects/omics_share/human/GRCh38/genome/annotation/function/COSMICv95_Coding_Noncoding.gatk_formatted.vcf' | Human Only - Used in variant annotation.

There are additional parameters that point to required scripts for processing data and computing metrics, which should not be changed.  

--filter_trim = '${projectDir}/bin/shared/filter_trim.py'  
--read_group_pyfile = '${projectDir}/bin/shared/read_group_from_fastq.py'  
--stats_agg = '${projectDir}/bin/wes/aggregate_stats_wes.py'  
--cosmic_annot = '${projectDir}/bin/shared/Cosmic_Annotation_hg38.pl'  





'''
}

