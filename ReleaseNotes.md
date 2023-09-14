# RELEASE NOTES

## Release 0.4.4

In this minor release we have corrected a syntax error in the parsing of single end CSV input to EMASE and GBRS. The syntax error prevented the workflow from running single end data when CSV input files were used.  

### Pipelines Added:

None

### Modules Added:

None

### Pipeline Changes:

1. EMASE: Correct csv single end parsing syntax.   
2. GBRS: Correct csv single end parsing syntax.  

### Module Changes:

None

## Release 0.4.3

In this minor release we have patched PTA to correct for a potential script error relating annotating CNVs and SVs on chromosome Y.  

### Pipelines Added:

None

### Modules Added:

None

### Pipeline Changes:

1. PTA: Adjusted when chromosome Y is included vs. excluded in caller merge and annotation steps.   

### Module Changes:

None

## Release 0.4.2

In this minor release we have made minor adjustments to the amplicon workflow, and added strandedness log output. 

### Pipelines Added:

None

### Modules Added:

None

### Pipeline Changes:

1. Amplicon: Alignment statistics are now taken post BQSR re-alignment.  

### Module Changes:

1. Primerclip: memory request increase.
2. python/python_check_strandedness.nf: added log file output.


## Release 0.4.1

In this release we have added one additional pipeline: amplicon sequencing. This pipeline support the analysis of [IDT xGen Amplicon panels](https://www.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-amplicon-sequencing/predesigned-amplicon-panels/), with current file support for [xGen Human Sample ID Amplicon Panel](https://www.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-amplicon-sequencing/predesigned-amplicon-panels/sample-id-amp-panel). Additionally, we have added a [classifier for EBV-associated PDX lymphomas](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6604205/) to the PDX RNA pipeline.

### Pipelines Added:

1. Amplicon

### Modules Added:

1. python/python_generate_fingerprint_report.nf  
2. python/python_lymphoma_classifier.nf  

### Pipeline Changes:

1. PDX RNAseq: added a classifier for EBV-associated PDX lymphomas.  

### Module Changes:

1. Cutadapt module function renamed from 'FILTER_FASTQ' to 'CUTADAPT'. Module file name adjusted to cutdadapt/cutadapt.nf   
2. python/python_check_strandedness.nf: Added strandedness override parameter for cases when `check_strandedness` fails to determine strand directionality. Corrected logic bug associated with parsing output from the tool.  
3. rsem/rsem_alignment_expression.nf: Resource request adjustment.  


## Release 0.4.0

In this release we have added five additional pipelines as part of the genetic diversity analysis suite. These pipelines support the analysis of genetically diverse samples (e.g., DO and CC mice) with [EMASE](https://github.com/churchill-lab/emase) and [GBRS](https://github.com/churchill-lab/GBRS), and the generation of reference files required for running these tools. 

### Pipelines Added:

1. EMASE
2. GBRS
3. Generate Pseudoreference
4. Prepare EMASE Reference/Inputs
5. Prepare DO GBRS Inputs

### Modules Added:

1. alntools/alntools_bam2emase.nf
2. bowtie/bowtie.nf
3. bowtie/bowtie_build.nf
4. emase/emase_create_hybrid.nf
5. emase/emase_get_common_alignment.nf
6. emase/emase_prepare_emase.nf
7. emase/emase_run.nf
8. g2gtools/g2gtools_convert.nf
9. g2gtools/g2gtools_extract.nf
10. g2gtools/g2gtools_gtf2db.nf
11. g2gtools/g2gtools_patch.nf
12. g2gtools/g2gtools_transform.nf
13. g2gtools/g2gtools_vcf2vci.nf
14. gbrs/gbrs_bam2emase.nf
15. gbrs/gbrs_compress.nf
16. gbrs/gbrs_export.nf
17. gbrs/gbrs_interpolate.nf
18. gbrs/gbrs_plot.nf
19. gbrs/gbrs_quantify.nf
20. gbrs/gbrs_quantify_genotype.nf
21. gbrs/gbrs_reconstruct.nf
22. python/append_dropped_chroms.nf
23. python/clean_prepEmase_transcriptList.nf
24. python/parse_gene_positions.nf
25. python/parse_transprobs.nf
26. r/do_transition_probablities.nf
27. r/generate_grid_file.nf
28. samtools/samtools_faidx_g2gtool.nf
29. utility_modules/filter_gtf_biotypes.nf
30. utility_modules/snorlax.nf

### Pipeline Changes:

None

### Module Changes:

None

## Release 0.3.1

In this minor release we have modified the behavior of Xenome to output compressed FASTQ files, and to delete the intermediate FASTQ files that are generated. We are implementing this change because the previous behavior of Xenome resulted in a large amount of redundant data in work directories.

We also added PDX test data for RNA-fusion.

### Pipelines Added:

None

### Modules Added:

None

### Pipeline Changes:

1. Changes to PDX RNA-seq, PDX WES, PDX RNA Fusion, and PDX PTA to reflect modifications to Xenome

### Module Changes:

1. xenome/xenome.nf modified to combine `xenome classify` and `fastq-sort` into the XENOME_CLASSIFY module. For non-fusion applications, human and mouse reads are now emitted as compressed .fastq.gz files
2. Removed fastq-tools/fastq-sort.nf as its functionality is now in xenome/xenome.nf
3. Modified input type specification for kallisto/kallisto_insert_size.nf to address issue with flash storage mounting in Singularity.
4. Added text file to pubDir statement in Picard collectRNAseqMetrics

## Release 0.3.0

In this major release we have added two additional pipelines, added flexibility for specifying inputs via sample sheets, support for downloading remote input data, support for GRCm39, support for PDX data, and many more changes detailed below. Additionally, we have added the concept of "subworkflows" for tasks that are more complex than a module and/or involve multiple containers, yet can be potentially re-used in multiple pipelines.

### Pipelines Added:

1. ChIP-seq - human, mouse
2. Paired Tumor Analysis (somatic/germline WGS) - human, PDX

### Subworkflows Added:

1. Aria download for remote input data
2. Concatenate paired tumor/normal FASTQ files
3. RNA-seq for PDX input data

### Modules Added:

1. arriba/arriba.nf
2. bamtools/bamtools_filter.nf
3. bcftools/bcftools_germline_filter.nf
4. bcftools/bcftools_intersect_lancet_candidates.nf
5. bcftools/bcftools_merge_callers.nf
6. bcftools/bcftools_remove_spanning.nf
7. bcftools/bcftools_split_multiallelic_regions.nf
8. bcftools/bcftools_split_multiallelic.nf
9. bedtools/bedtools_amplicon_metrics.nf
10. bedtools/bedtools_genomecov.nf
11. bedtools/bedtools_start_candidates.nf
12. biqseq2/bicseq2_normalize.nf
13. biqseq2/bicseq2_seg_unpaired.nf
14. biqseq2/bicseq2_seg.nf
15. conpair/conpair_pileup.nf
16. conpair/conpair.nf
17. cosmic/cosmic_add_cancer_resistance_mutations_germline.nf
18. cosmic/cosmic_add_cancer_resistance_mutations_somatic.nf
19. cosmic/cosmic_annotation_somatic.nf
20. cosmic/cosmic_annotation.nf
21. deeptools/deeptools_computematrix.nf
22. deeptools/deeptools_plotfingerprint.nf
23. deeptools/deeptools_plotheatmap.nf
24. deeptools/deeptools_plotprofile.nf
25. ensembl/varianteffectpredictor_germline.nf
26. ensembl/varianteffectpredictor_somatic.nf
27. fastq-tools/fastq-pair.nf
28. fastq-tools/fastq-sort.nf
29. fusion_report/fusion_report.nf
30. fusioncatcher/fusioncatcher.nf
31. gatk/gatk_cnnscorevariants.nf
32. gatk/gatk_combinegvcfs.nf
33. gatk/gatk_filtermutectcalls_tumorOnly.nf
34. gatk/gatk_filtermutectcalls.nf
35. gatk/gatk_filtervarianttranches.nf
36. gatk/gatk_genotype_gvcf.nf
37. gatk/gatk_getsamplename_noMeta.nf
38. gatk/gatk_getsamplename.nf
39. gatk/gatk_haplotypecaller_sv_germline.nf
40. gatk/gatk_mergemutectstats.nf
41. gatk/gatk_mutect2_tumorOnly.nf
42. gatk/gatk_mutect2.nf
43. gatk/gatk_sortvcf_germline.nf
44. gatk/gatk_sortvcf_somatic_merge.nf
45. gatk/gatk_sortvcf_somatic_tools.nf
46. gatk/gatk_variantfiltration_af.nf
47. gatk/gatk_variantfiltration_mutect2.nf
48. gatk/gatk3_applyrecalibration.nf
49. gatk/gatk3_genotypegvcf.nf
50. gatk/gatk3_haplotypecaller.nf
51. gatk/gatk3_indelrealigner.nf
52. gatk/gatk3_realignertargetcreator.nf
53. gatk/gatk3_variantannotator.nf
54. gatk/gatk3_variantrecalibrator.nf
55. gridss/gridss_assemble.nf
56. gridss/gridss_calling.nf
57. gridss/gridss_chrom_filter.nf
58. gridss/gridss_preprocess.nf
59. gridss/gripss_somatic_filter.nf
60. homer/annotate_boolean_peaks.nf
61. homer/homer_annotatepeaks.nf
62. homer/plot_homer_annotatepeaks.nf
63. illumina/manta.nf
64. illumina/strelka2.nf
65. jaffa/jaffa.nf
66. kallisto/kallisto_insert_size.nf
67. kallisto/kallisto_quant.nf
68. lumpy_sv/lumpy_sv.nf
69. macs2/macs2_consensus.nf
70. macs2/macs2_peak_calling_chipseq.nf
71. macs2/plot_macs2_qc.nf
72. msisensor2/msisensor2_tumorOnly.nf
73. msisensor2/msisensor2.nf
74. multiqc/multiqc_custom_phantompeakqualtools.nf
75. novocraft/novosort.nf
76. nygc-short-alignment-marking/short_alignment_marking.nf
77. nygenome/lancet_confirm.nf
78. nygenome/lancet.nf
79. phantompeakqualtools/phantompeakqualtools.nf
80. picard/picard_cleansam.nf
81. picard/picard_collectmultiplemetrics.nf
82. picard/picard_collecttargetpcrmetrics.nf
83. picard/picard_fix_mate_information.nf
84. picard/picard_mergesamfiles.nf
85. pizzly/pizzly.nf
86. preseq/preseq.nf
87. primerclip/primerclip.nf
88. python/python_add_final_allele_counts.nf
89. python/python_add_nygc_allele_counts.nf
90. python/python_check_strandedness.nf
91. python/python_filter_pon.nf
92. python/python_filter_vcf.nf
93. python/python_germline_vcf_finalization.nf
94. python/python_get_candidates.nf
95. python/python_merge_columns.nf
96. python/python_merge_prep.nf
97. python/python_remove_contig.nf
98. python/python_rename_metadata.nf
99. python/python_rename_vcf.nf
100. python/python_reorder_vcf_columns.nf
101. python/python_snv_to_mnv_final_filter.nf
102. python/python_somatic_vcf_finalization.nf
103. python/python_split_mnv.nf
104. python/python_vcf_to_bed.nf
105. r/annotate_bicseq2_cnv.nf
106. r/annotate_genes_sv.nf
107. r/annotate_sv_with_cnv.nf
108. r/annotate_sv.nf
109. r/filter_bedpe.nf
110. r/frag_len_plot.nf
111. r/merge_sv.nf
112. samtools/samtools_faidx.nf
113. samtools/samtools_filter_unique_reads.nf
114. samtools/samtools_filter.nf
115. samtools/samtools_mergebam_filter.nf
116. samtools/samtools_stats_insertsize.nf
117. samtools/samtools_stats.nf
118. samtools/samtools_view.nf
119. squid/squid_annotate.nf
120. squid/squid_call.nf
121. star/star_align.nf
122. star-fusion/star-fusion.nf
123. subread/subread_feature_counts_chipseq.nf
124. svaba/svaba.nf
125. tabix/compress_merged_vcf.nf
126. tabix/compress_vcf_region.nf
127. tabix/compress_vcf.nf
128. ucsc/ucsc_bedgraphtobigwig.nf
129. utility_modules/aria_download.nf
130. utility_modules/chipseq_bampe_rm_orphan.nf
131. utility_modules/chipseq_check_design.nf
132. utility_modules/chipseq_make_genome_filter.nf
133. utility_modules/concatenate_reads_sampleSheet.nf
134. utility_modules/deseq2_qc.nf
135. utility_modules/frip_score.nf
136. utility_modules/get_read_length.nf
137. utility_modules/gunzip.nf
138. utility_modules/jax_trimmer.nf
139. utility_modules/parse_extracted_sv_table.nf
140. xenome/xenome.nf

### Pipeline Changes:

1. WES, RNA-seq, and RNA-fusion added support for PDX data  
2. WES, RNA-seq, WGS, ATAC, RRBS, ChIP added support for GRCm39  
3. Support for input specification using sample sheets for ATAC, RNA-seq, RRBS, WES, WGS  
4. Support for downloading input data for ATAC, RNA-seq, RRBS, WES, WGS  
5. Added MULTIQC to ATAC, RNA-seq, RRBS, WES, WGS  
6. Added assessment of strandedness using python/python_check_strandedness.nf rather than requiring specification via parameters  
7. Added assessment of read length for RNAseq for STAR index selection rather than requiring specfication via parameters  
8. Modified variant annotations in WES and WGS  
9. Added GVCF support for WES and WGS  

### Module Changes:

1. errorStrategy modified for all modules to catch and report instances where tasks fail due to walltime or memory contraints. This previously required a deep reading of the subtask SLURM logs, but now will be reported in the top-level SLURM log and is more user-friendly
2. Removed log.info statements from modules to avoid noisy disruption of log files
3. ChIP-seq support for bwa/bwa_mem.nf, fastqc/fastqc.nf, picard/picard_markduplicates.nf, trim_galore/trim_galore.nf
4. Corrected emit statements for g2gtools/g2gtools_chain_convert_peak.nf
5. Corrected emit statements for gatk/gatk_chain_filter_reads.nf
6. Modified gatk/gatk_haplotypecaller_interval.nf and gatk/gatk_haplotypecaller.nf for optional GVCF support
7. Generalized multiqc/multiqc.nf via parameter for multiqc config
8. Removed --METRIC_ACCUMULATION_LEVEL ALL_READS and --VALIDATION_STRINGENCY LENIENT parameters from picard/picard_collectalignmentsummarymetrics.nf
9. Modified strand specification logic for picard/picard_collectrnaseqmetrics.nf
10. Updated rsem/rsem_alignment_expression.nf to reflect changes in strandedness detection, reorganized outputs and catching log files for multiqc
11. Changes to output text for mt DNA content in samtools/samtools_calc_mtdna_filter_chrm.nf
12. Changes to output text from samtools/samtools_final_calc_frip.nf
13. Changes to output formatting for samtools/samtools_quality_checks.nf
14. Updated snpEff container to v5.1d to support GRCm39
15. Changes to output fields for mouse and human from snpeff_snpsift/snpsift_extractfields.nf
16. Added missing container to utility_modules/concatenate_reads_PE.nf and utility_modules/concatenate_reads_SE.nf

## Release 0.2.2

* Change WES and WGS COMSIC annotation to use SNPsift. 
* Added explicit dbSNP annotation. 

### Pipelines Added:

NONE

### Modules Added:

1. SNPSIFT_ANNOTATE 

### Pipeline Changes:

1. WES and WGS now use SNPSift to annotate COSMIC and dbSNP IDs onto variants. 

### Module Changes:

1. COSMIC_ANNOTATION and associated perl scripts removed. 


## Release 0.2.1

Added STAR support to RNA-seq pipeline.

### Pipelines Added:

NONE

### Modules Added:

NONE

### Pipeline Changes:

1. RNA-seq pipeline now supports STAR and bowtie2 (default) through the RSEM module.

### Module Changes:

1. RSEM: --rsem_aligner accepts "bowtie2" or "star." The default STAR indices for mouse and human are 100 bp, with alternates suggested in the RNA-seq config file.

## Release 0.2.0

**NOTE:** This release contains a patch for multi-sample processing. We strongly recommend multi-sample processing done prior to this release should be re-run with v0.2.0+

### Pipelines Added:

1. RRBS - Mouse & Human
2. ATAC - Mouse & Human

### Modules Added: 

1. FastQC
2. Trim-Galore
3. Bismark Alignment
4. Bismark Deduplicator
5. Bismark Methylation Extractor 
6. MultiQC
7. Bedtools functions for ATAC QC summary
8. Bowtie2
9. Cutadapt
10. Deeptools bamcoverage and alignmentSieve
11. g2gTools chain convert
12. Macs2 ATAC peak calling and ATAC peak coverage
13. Subread feature counts

### Pipeline Changes:

1. Multiple pipeline changes related to multi-sample patch.
2. Modified module load statements to invoke "${projectDir}" instead of relative "../" path.
3. Removed CTP and Probe coverage calculations from human RNA-seq

### Module Changes:

1. Multiple module changes related to multi-sample patch. 
2. Trimmomatic Trim stub module removed. 
3. RSEM - forward stranded option added. 
4. Picard Collect RNAseqMetrics - forward strand option added. 

## Release 0.1.2 

Updated run scripts to load CS supported Nextflow module. 

## Release 0.1.1 

### Pipelines Added:

NONE

### Modules Added: 

1. concatenate_reads_PE.nf
2. concatenate_reads_SE.nf
3. Modules refactored to individual files (e.g., gatk_haplotypecaller.nf). 

### Pipeline Changes:

1. Added ability to concatenate Fastq files by sample, which are split across sequencing lanes into single R1/R2 or R1 files (depending on PE or SE). 
2. Adjusted pipelines for refactored module files.
3. Fixed CTP/PROBE typo in human RNA coverage calculation.
4. Added HPC `--profile` options and settings for Sumner and Elion. 

### Module Changes:

1. Adjusted WGS wall clock settings. 
1. Refactored modules to individual files (e.g., gatk_haplotypecaller.nf). 
2. Set pipeline script parameter to hard coded paths.
3. Cleaned all Nextflow files from the bin directory.
4. Removed Sumner specific HPC settings from each module. 

## Release 0.1.0 -- 03.28.2022

### Pipelines Added:
1. Whole Genome Sequencing - Mouse & Human
2. Whole Exome Sequencing - Mouse & Human
3. RNA Sequencing - Mouse & Human

### Modules Added: 
1. bamtools.nf
1. bcftools.nf
1. bwa.nf
1. cosmic.nf
1. gatk.nf
1. picard.nf
1. quality_stats.nf
1. read_groups.nf
1. rsem.nf
1. samtools.nf
1. snpeff.nf
1. snpsift.nf
1. summary_stats.nf
1. trimmomatic.nf

### Pipeline Changes:
NONE

### Module Changes:
NONE