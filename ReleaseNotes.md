# RELEASE NOTES

## Release 0.6.7

In this release we make the following minor adjustments:   

1. Correct syntax errors in the Xengsort module when running single-end data.  
1. Minor adjustments to EMASE and GBRS help and log information to include the `gen_org` param.  
1. Bump the version of MultiQC to v1.23.  
1. Increase the memory request for a `PTA` moudles: `python_merge_prep.nf` and `python_reorder_vcf_columns.nf`.  
1. Add `CHECK_STRANDEDNESS` to multiQC output for PDX RNAseq
1. Increased job memory request in example run scripts.   

## Release 0.6.6

In this release, we add a FASTQ sorting function to the Xengsort module. Due to asynchronous multi-threading in the classification step, Xengsort produces FASTQ output with non-deterministic sort order. BWA produces subtly different mapping results when reads in otherwise identical FASTQ inputs are shuffled ([see note from BWA developer here](https://github.com/lh3/bwa/issues/192#issuecomment-380612006)). The slight mapping differences are not enough to impact overall results, but do prevent fully reproducible results when Xengsort is used and reads are not sorted. The addition of the sorting function allows for fully reproducible results, with no additional user action required.    

## Release 0.6.5

In this minor release, we fix a `subscript out of bounds` bug in `bin/wes/sequenza_seg_na_window.R`.    

## Release 0.6.4

In this release, we adjust memory and wallclock requirements for a number of modules, update `read_group_from_fastq.py` from python2 to python3, and incorporate PRs #4 and #5.  

* PR #4 (contributed by @BrianSanderson) adds an optional gene and transcript count merge across samples in the RNA and PDX RNA workflows (merge accessed via including the `--merge_rna_counts` flag).  
* PR #5 (contributed by @alanhoyle) adds a catch for corrupt gzip files in the Bowtie module as used by EMASE/GRBS analyses.  

### Pipelines Added:

None

### Modules Added:

1. utility_modules/merge_rsem_counts.nf

### Pipeline Changes:

1. workflows/rnaseq.nf module added to merge gene and transcript expression when `--merge_rna_counts` is used.  
1. workflows/pdx_rnaseq.nf module added to merge gene and transcript expression when `--merge_rna_counts` is used.  

### Module Changes:

1. bowtie/bowtie.nf pipefail catch added for corrupt gzip files, per #5. 
1. fastp/fastp.nf save json report as well as html report.  
1. nygenome/lancet.nf wallclock request increase.  
1. picard/picard_markduplicates.nf memory adjustment, and accounting for MarkDuplicates not fully respecting -Xmx memory limits imposed by Java.  
1. picard/picard_reordersam.nf memory request increase.  
1. picard/picard_sortsam.nf memory request increase.  
1. utility_modules/read_groups.nf container changed to py3.  

### Script Changes:

1. bin/shared/read_group_from_fastq.py update from py2 to py3. 


## Release 0.6.3

In this release we change the read disambiguation tool Xenome for Xengsort. Extensive benchmarking shows high concordance among results obtained from both tools.  

Additionally, we correct an issue with the human PTA workflow when running the combination of the `--pdx` and `--split_fastq` options. Data run with this combination of options from version 0.6.0-0.6.2 should be re-run. 

### Pipelines Added:

None

### Modules Added:

1. xengsort/xengsort_classify.nf
1. xengsort/xengsort_index.nf

### Pipeline Changes:

1. Xengsort replaces Xenome for all PDX based workflows (RNAseq, RNA fusion, Hs PTA, Somatic WES, Somatic WES PTA)
1. Correction made for the Human PTA when running the combination of the `--pdx` and `--split_fastq` options.

### Module Changes:

None


## Release 0.6.2

In this minor release we adjust memory and wall clock statements, and modified `bin/pta/merge-caller-vcfs.r` to correct for an edge case related bug.

## Release 0.6.1

In this minor release we added support for automatic Zenodo releases via github actions. There are no changes or additions to workflows.

## Release 0.6.0

In this major release we add seven new workflows, and make numerous changes to existing workflows. Specific changes are discussed below.  

For Jackson Laboratory users this release now supports the new Sumner2 cluster. To use workflows on sumner, simply specify: `-profile sumner2`. Note that Sumner has reached end of life, and will no longer be supported going forward. We have updated all example run scripts to use `-profile sumner2`.  

**Note** Sumner2 enforces strict Linux cgroups, which holds jobs to the memory and cpu limits requested by each Nextflow module. In our release testing, we increased many memory reservation steps; however, additional memory issues are to be expected. If you encounter `OOM` (out of memory) issues and experience workflow steps failing with `killed` reported in the error log, please either email us: (ngsOps@jax.org) or submit an [issue](https://github.com/TheJacksonLaboratory/cs-nf-pipelines/issues) with details on which module failed and the size of the dataset you were running. 

Related to memory and time restrictions, we made signficiant changes to the PTA, WGS, and WES workflows:    

1. For human PTA, WGS and WES analyses GATK BaseRecalibrator is now scattered by chromosome.   
1. For PTA and WGS, options were added to allow users to:  
    1. Deduplicate reads with [`Clumpify`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/clumpify-guide/) prior to mapping steps. 
    1. Split FASTQ files into batched chunks for subsequent mapping. Mapped batches are merged prior to the GATK MarkDuplicates step.  
    1. Cap coverage at a user defined threshold using [JVARKIT Biostar154220](http://lindenb.github.io/jvarkit/Biostar154220.html) prior to variant calling. This can help reduce computational load when calling variants in higher coverage areas of the genome. 
1. For PTA, WGS and WES FASTP is now used for read and adapter trimming. 

We have included an option to specify the root location of the omics_share reference file set. For Jax users on Sumner2, this option should not be changed and defaults to `/projects/omics_share`. For external users, or those on Elion, specify the root directory of reference files with `--reference_cache </path/to/omics_share>`.   

Finally, we added testing modules for all workflows to be used with (nf-test)[https://www.nf-test.com/].  

### Pipelines Added:

1. Amplicon Sequencing (supporting human only at this time): General PCR / Targeted Sequencing
1. Genetic Ancestry (See https://www.biorxiv.org/content/10.1101/2022.10.24.513591v1 for details and methods)
1. Germline Structural Variant Calling
    1. Illumina short-read data
    1. Pacific Biosciences (PacBio) long-read data: CCS and CLR modes
    1. Oxford Nanopore Technologies long-read data (ONT)
1. Somatic Whole Exome Sequencing for tumor-only samples (with option for PDX)
1. Somatic Whole Exome Sequencing for Paired Tumor Analysis (PTA; with option for PDX)

### Modules Added:

1. modules/abra2/abra2.nf
1. modules/bbmap/bbmap_clumpify.nf
1. modules/bcftools/bcftools_annotate.nf
1. modules/bcftools/bcftools_call.nf
1. modules/bcftools/bcftools_duphold_filter.nf
1. modules/bcftools/bcftools_filter.nf
1. modules/bcftools/bcftools_merge_amplicon.nf
1. modules/bcftools/bcftools_mpileup.nf
1. modules/bcftools/bcftools_norm.nf
1. modules/bcftools/bcftools_rehead_sort.nf
1. modules/bcftools/bcftools_vcf_to_bcf.nf
1. modules/bedops/bedops_sort.nf
1. modules/bedops/bedops_window.nf
1. modules/bedtools/bedtools_sequenza_subtract.nf
1. modules/bwa/bwa_index.nf
1. modules/bwa/bwa_mem2.nf
1. modules/delly/delly_call.nf
1. modules/delly/delly_call_germline.nf
1. modules/delly/delly_cnv_germline.nf
1. modules/duphold/duphold.nf
1. modules/freebayes/freebayes.nf
1. modules/gatk/gatk_baserecalibrator_interval.nf
1. modules/gatk/gatk_calculatecontamination.nf
1. modules/gatk/gatk_calculatecontamination_tumorOnly.nf
1. modules/gatk/gatk_filtermutectcalls_wes.nf
1. modules/gatk/gatk_gatherbqsrreports.nf
1. modules/gatk/gatk_getpileupsummaries.nf
1. modules/gatk/gatk_getpileupsummaries_tumorOnly.nf
1. modules/gatk/gatk_haplotypecaller_amplicon.nf
1. modules/gatk/gatk_learnreadorientationmodel.nf
1. modules/gatk/gatk_mutect2_wes_pta.nf
1. modules/gatk/gatk_printreads.nf
1. modules/gatk/gatk_variantfiltration_freebayes.nf
1. modules/illumina/manta_germline.nf
1. modules/jvarkit/jvarkit_biostar154220.nf
1. modules/lumpy/lumpy_call_sv.nf
1. modules/lumpy/lumpy_extract_splits.nf
1. modules/lumpy/lumpy_prep.nf
1. modules/minimap/minimap2_index.nf
1. modules/minimap/minimap2_map_ont.nf
1. modules/nanofilt/nanofilt.nf
1. modules/nanoqc/nanoqc.nf
1. modules/nanostat/nanostat.nf
1. modules/nanosv/nanosv.nf
1. modules/pbmm2/pbmm2_call.nf
1. modules/pbmm2/pbmm2_index.nf
1. modules/pbsv/pbsv_call.nf
1. modules/pbsv/pbsv_discover.nf
1. modules/picard/picard_markduplicates_removedup.nf
1. modules/picard/picard_sortsam_mmrsvd.nf
1. modules/porechop/porechop.nf
1. modules/python/python_add_AF_freebayes.nf
1. modules/python/python_add_AF_haplotypecaller.nf
1. modules/python/python_annot_depths.nf
1. modules/python/python_annot_on_target.nf
1. modules/python/python_bedpe_to_vcf.nf
1. modules/python/python_parse_depths.nf
1. modules/python/python_parse_survivor_ids.nf
1. modules/r/illumina_sv_merge.nf
1. modules/r/r_merge_depths.nf
1. modules/samtools/samtools_cat.nf
1. modules/samtools/samtools_filter_mmrsvd.nf
1. modules/samtools/samtools_merge.nf
1. modules/samtools/samtools_mpileup.nf
1. modules/samtools/samtools_stats_mmrsvd.nf
1. modules/scarhrd/scarhrd.nf
1. modules/sequenza/sequenza_annotate.nf
1. modules/sequenza/sequenza_na_window.nf
1. modules/sequenza/sequenza_pileup2seqz.nf
1. modules/sequenza/sequenza_run.nf
1. modules/smoove/smoove_call_germline.nf
1. modules/sniffles/sniffles.nf
1. modules/snpweights/snpweights_inferanc.nf
1. modules/snpweights/snpweights_vcf2eigenstrat.nf
1. modules/survivor/survivor_annotation.nf
1. modules/survivor/survivor_bed_intersect.nf
1. modules/survivor/survivor_inexon.nf
1. modules/survivor/survivor_merge.nf
1. modules/survivor/survivor_to_bed.nf
1. modules/survivor/survivor_vcf_to_table.nf
1. modules/tumor_mutation_burden/tmb_score.nf
1. modules/utility_modules/filter_trim.nf
1. modules/vcftools/vcftools_filter.nf

### NF-Test Modules Added: 

1. tests/workflows/amplicon_fingerprint.nf.test
1. tests/workflows/amplicon_generic.nf.test
1. tests/workflows/ancestry.nf.test
1. tests/workflows/atac.nf.test
1. tests/workflows/chipseq.nf.test
1. tests/workflows/emase.nf.test
1. tests/workflows/gbrs.nf.test
1. tests/workflows/generate_pseudoreference.nf.test
1. tests/workflows/prep_do_gbrs_inputs.nf.test
1. tests/workflows/prepare_emase.nf.test
1. tests/workflows/pta.nf.test
1. tests/workflows/rna_fusion.nf.test
1. tests/workflows/rnaseq.nf.test
1. tests/workflows/rrbs.nf.test
1. tests/workflows/somatic_wes.nf.test
1. tests/workflows/somatic_wes_pta.nf.test
1. tests/workflows/wes.nf.test
1. tests/workflows/wgs.nf.test

### Pipeline Changes:

1. chipseq.nf: Error reporting added for malformed CSV input files  
1. pta.nf: Error reporting added for malformed CSV input files  
1. subworkflows/hs_pta.nf: `JAX_TRIMMER` replaced with `FASTP`. GATK Baserecalibration is now scattered by chromosome. Options added to: 1. deduplicate reads with [`Clumpify`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/clumpify-guide/) prior to mapping steps, 2. Split FASTQ files into batched chunks for subsequent mapping. 3. Cap coverage at a user defined threshold using [JVARKIT Biostar154220](http://lindenb.github.io/jvarkit/Biostar154220.html) prior to variant calling. Additionally, short_alignment_marking following mapping was previously disconnected for the workflow. This step has been included.  
1. subworkflows/mm_pta.nf: `JAX_TRIMMER` replaced with `FASTP`. Options added to: 1. deduplicate reads with [`Clumpify`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/clumpify-guide/) prior to mapping steps, 2. Split FASTQ files into batched chunks for subsequent mapping. 3. Cap coverage at a user defined threshold using [JVARKIT Biostar154220](http://lindenb.github.io/jvarkit/Biostar154220.html) prior to variant calling.
1. rnaseq.nf: `Check Strandedness` log data added to MultiQC report.  
1. wes.nf: `JAX_TRIMMER` replaced with `FASTP`. For human analysis, GATK Baserecalibration is now scattered by chromosome.  
1. wgs.nf: `JAX_TRIMMER` replaced with `FASTP`. For human analysis, GATK Baserecalibration is now scattered by chromosome. Options added to: 1. deduplicate reads with [`Clumpify`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/clumpify-guide/) prior to mapping steps, 2. Split FASTQ files into batched chunks for subsequent mapping. 3. Cap coverage at a user defined threshold using [JVARKIT Biostar154220](http://lindenb.github.io/jvarkit/Biostar154220.html) prior to variant calling.  


### Module Changes:

1. alntools/alntools_bam2emase.nf: Bump docker container version.  
1. bedtools/bedtools_genomecov.nf: Memory request increase.  
1. biqseq2/bicseq2_normalize.nf: Adjustment to read length parsing logic.  
1. bowtie/bowtie.nf: Memory request increase.  
1. bwa/bwa_mem.nf: Input tuple adjustment.  
1. bwa/bwa_mem_hla.nf: Input tuple adjustment. 
1. deeptools/deeptools_filter_remove_multi_sieve.nf
1. emase/emase_create_hybrid.nf: Bump docker container version.  
1. emase/emase_get_common_alignment.nf: Bump docker container version. Memory request increase.  
1. emase/emase_prepare_emase.nf: Bump docker container version.  
1. emase/emase_run.nf: Bump docker container version.  
1. ensembl/varianteffectpredictor_germline_mouse.nf: Input tuple adjustment. Add BGZIP and indexing to final VCF output.  
1. fastp/fastp.nf: Memory request increase.  
1. gatk/gatk_applybqsr.nf: Added tmp dir to command.  
1. gatk/gatk_baserecalibrator.nf: Added tmp dir to command.  
1. gatk/gatk_chain_extract_badreads.nf: Added tmp dir to command.  
1. gatk/gatk_chain_filter_reads.nf: Added tmp dir to command.  
1. gatk/gatk_cnnscorevariants.nf: Added tmp dir to command.  
1. gatk/gatk_combinegvcfs.nf: Added tmp dir to command.  
1. gatk/gatk_depthofcoverage.nf: Added tmp dir to command.  
1. gatk/gatk_filtermutectcalls.nf: Added tmp dir to command.  
1. gatk/gatk_filtervarianttranches.nf: Added tmp dir to command.  
1. gatk/gatk_genotype_gvcf.nf: Added tmp dir to command.  
1. gatk/gatk_getsamplename.nf: Added tmp dir to command.  
1. gatk/gatk_getsamplename_noMeta.nf: Added tmp dir to command.  
1. gatk/gatk_haplotypecaller.nf: Added tmp dir to command.  
1. gatk/gatk_haplotypecaller_interval.nf: Added tmp dir to command.  
1. gatk/gatk_haplotypecaller_sv_germline.nf: Added tmp dir to command.  
1. gatk/gatk_indexfeaturefile.nf: Added tmp dir to command.  
1. gatk/gatk_mergemutectstats.nf: Added tmp dir to command.  
1. gatk/gatk_mergevcf.nf: Added tmp dir to command.  
1. gatk/gatk_mergevcf_list.nf: Added tmp dir to command.  
1. gatk/gatk_mutect2.nf: Added tmp dir to command.  
1. gatk/gatk_mutect2_tumorOnly.nf: Added tmp dir to command.  
1. gatk/gatk_selectvariants.nf: Added tmp dir to command.  
1. gatk/gatk_sortvcf_germline.nf: Added tmp dir to command.  
1. gatk/gatk_sortvcf_somatic_merge.nf: Added tmp dir to command.  
1. gatk/gatk_sortvcf_somatic_tools.nf: Added tmp dir to command.  
1. gatk/gatk_updatevcfsequencedictionary.nf: Added tmp dir to command.  
1. gatk/gatk_variantfiltration.nf: Added tmp dir to command.  
1. gatk/gatk_variantfiltration_af.nf: Added tmp dir to command.  
1. gatk/gatk_variantfiltration_mutect2.nf: Added tmp dir to command.  
1. gbrs/gbrs_bam2emase.nf: Bump docker container version. Memory request increase.  
1. gbrs/gbrs_compress.nf: Bump docker container version. Memory request increase.  
1. gbrs/gbrs_export.nf: Bump docker container version. 
1. gbrs/gbrs_interpolate.nf: Bump docker container version. 
1. gbrs/gbrs_plot.nf: Bump docker container version. 
1. gbrs/gbrs_quantify.nf: Bump docker container version. 
1. gbrs/gbrs_quantify_genotype.nf: Bump docker container version. 
1. gbrs/gbrs_reconstruct.nf: Bump docker container version. 
1. gridss/gridss_assemble.nf: Memory request increase.  
1. illumina/strelka2.nf: Wallclock request increase.  
1. multiqc/multiqc.nf: Tool version updated to v1.21
1. nygc-short-alignment-marking/short_alignment_marking.nf: Bug correction in original module script.  
1. picard/picard_cleansam.nf: Output naming adjustment. 
1. picard/picard_collectalignmentsummarymetrics.nf: Added tmp dir to command.  
1. picard/picard_collecttargetpcrmetrics.nf: Added tmp dir to command.  
1. picard/picard_collectwgsmetrics.nf: Added tmp dir to command.  
1. picard/picard_fix_mate_information.nf: Corrected BAM sort order of output to coordinate.  
1. picard/picard_sortsam.nf: Added index creation option.  
1. samtools/samtools_calc_mtdna_filter_chrm.nf: Memory request increase.  
1. samtools/samtools_faidx.nf: Input tuple adjustment, and output reorganization.  
1. snpeff_snpsift/snpeff_snpeff.nf: Memory request increase. Tmp dir adjustment.  
1. snpeff_snpsift/snpsift_extractfields.nf: Added support for amplicon_generic, somatic_wes, and somatic_wes_pta workflows.  
1. squid/squid_call.nf: Memory request increase.  
1. utility_modules/chipseq_make_genome_filter.nf: Input tuple adjustment.  
1. utility_modules/jax_trimmer.nf: File output naming adjusted.  

### Script Added: 

ancestry/vcf2eigenstrat.py: Convert VCF to EigenStrat format.  
germline_sv/annot_vcf_with_depths.py: Add info fields for depths from individual caller to VCF files.  
germline_sv/annot_vcf_with_exon.py: Apply 'InExon' INFO fields to original SV VCF files.  
germline_sv/annot_vcf_with_on_target.py: Apply 'OnTarget' tINFO fields to original SV VCF files.  
germline_sv/bedpetovcf.py: Convert BEDPE format back to `SURVIVOR` like VCF.  
germline_sv/clean_sniffles.sh: Adjust `Sniffles` calls.  
germline_sv/cnvnator2VCF.pl: Convert `CNVnator` formatted files to VCF.  
germline_sv/hydra_to_vcf.py: Convert Hydra BEDPE output into VCF 4.1 format.  
germline_sv/merge_depths.R: Merge nanoSV and Sniffles read/support depths.  
germline_sv/merge_sv.r: Merge an arbitrary number of VCFs, and annotate with simple event type.  
germline_sv/parse_caller_depths.py: Parse SV caller VCFs to extract IDs and depth information.  
germline_sv/parse_survivor_ids.py: Parse SURVIVOR merged VCFs to extract IDs.  
germline_sv/sed_unquote.sh: script to remove double-quotes from files, which is used avoid issues with unescaped quotes in Nextflow script blocks.  
germline_sv/summarize_intersections.R: Intersect SV calls by type with known structural variant databases.  
germline_sv/surv_annot.sh: Adjust `SURVIVOR` output to txt.   
germline_sv/surv_annot_process.R: Adjust surv_annot output by SV type.  
germline_sv/sv_to_table.py: Parse SURVIVOR merged VCF to output a summary table for each variant that lists the position, type, and size.  
wes/AF_freebayes.py: Add Estimated Allele Frequency (ALT_AF) to the INFO field of FreeBayes VCF output.  
wes/AF_haplotypecaller.py: Add Estimated Allele Frequency (ALT_AF) to the INFO field of HaplotypeCaller VCF output.  
wes/TMB_calc.R: Compute tumor mutation burden. See Somatic WES wiki for details.  
wes/allele_depth_min_and_AF_from_ADs.py: ecompute the locus depth from the allele-depths, and filter based on a minimum total allele depth. 
wes/ensembl_annotation.pl: Annotates Ensembl transcripts and genes with copy number and breakpoints.  
wes/scarHRD.R: Compute homologous recombination deficiency (HRD) with scarHRD.  
wes/sequenza_run.R: Compute copy number variantion with Sequenza.  
wes/sequenza_seg_na_window.R: Filter Sequenza CNV segments with `NA` calls within 1Mb windows. 

### Script Changes: 

gbrs/gene_bp_to_cM_to_transprob.R: Added local BIOMART_CACHE location.  
pta/make_main_vcf.py: Adjusted genomic build check logic blocks.  
pta/make_txt.py: Adjusted genomic build check logic blocks.  
pta/merge-caller-vcfs.r: Added logic to catch edge case where no variants were within a VCF for merging.  
shared/extract_csv.nf: Added error reporting for malformed CSV input files.  
shared/extract_gbrs_csv.nf: Added error reporting for malformed CSV input files.  

## Release 0.5.0

In this release we have added the mouse version of PTA, and changed the read trimmer for the RNAseq pipeline to Fastp. Additionally, the latest version of Nextflow is now supported.

Note for Jackson Laboratory users on the Sumner cluster: Fastscratch has reached end of life, and is no longer supported. We have updated all example run scripts to point at `/flashscratch` rather than `/fastscratch`. For production analyses all working directories (i.e., `-w <PATH>`) should use `/flashscratch/$USER/...`. 

### Pipelines Added:

1. Mouse PTA

### Modules Added:

1. bcftools/bcftools_bcf_to_vcf.nf
2. bcftools/bcftools_compress_index.nf
3. bcftools/bcftools_merge_delly_cnv.nf
4. bcftools/bcftools_query_delly_cnv.nf
5. delly/delly_call_somatic.nf
6. delly/delly_classify.nf
7. delly/delly_cnv_somatic.nf
8. delly/delly_filter_somatic.nf
9. ensembl/varianteffectpredictor_germline_mouse.nf
10. ensembl/varianteffectpredictor_somatic_mouse.nf
11. fastp/fastp.nf
12. gatk/gatk_updatevcfsequencedictionary.nf
13. python/python_somatic_vcf_finalization_mouse.nf
14. r/annotate_delly_cnv.nf
15. r/annotate_genes_sv_mouse.nf
16. r/annotate_sv_mouse.nf
17. r/annotate_sv_with_cnv_mouse.nf
18. r/filter_bedpe_mouse.nf
19. r/merge_sv_mouse.nf
20. r/plot_delly_cnv.nf
21. smoove/smoove_call.nf
22. svtyper/svtyper.nf
23. utility_modules/gzip.nf
24. utility_modules/lumpy_compress_index.nf

### Pipeline Changes:

1. RNAseq: The read trimmer script was replaced with `fastp`. STAR logs from RSEM now saved and passed to MultiQC for summary.
2. Human PTA: The read trimmer script was replace with `fastp`.

### Module Changes:

1. bwa/bwa_mem.nf: Wallclock and memory request adjustment.
2. emase/emase_get_common_alignment.nf: Wallclock request adjustment.
3. gatk/gatk_applybqsr.nf: Wallclock request adjustmnet.
4. gatk/gatk_sortvcf_somatic_tools.nf: Added mouse PTA support.
5. gridss/gridss_assemble.nf: Update container to correct bug in prior container build. Wallclock and memory adjustment.
6. gridss/gridss_calling.nf: Update container to correct bug in prior container build.
7. gridss/gridss_preprocess.nf: Update container to correct bug in prior container build.
8. lumpy_sv/lumpy_sv.nf: Modified previously unused module for use in mouse PTA.
9. msisensor2/msisensor2.nf: Correct `cp` error that can occur on nextflow resume.
10. msisensor2/msisensor2_tumorOnly.nf: Correct `cp` error that can occur on nextflow resume.
11. multiqc/multiqc.nf: Added cpu, memory, and wallclock requests.
12. nygenome/lancet.nf: Memory request adjustment.
13. nygenome/lancet_confirm.nf: Memory request adjustment.
14. picard/picard_addorreplacereadgroups.nf: Memory request adjustment. Adjusted PICARD temp directory to Nextflow work directory.
15. picard/picard_collectalignmentsummarymetrics.nf: Wallclock request adjustment.
16. picard/picard_collecthsmetrics.nf: Wallclock request adjustment.
17. picard/picard_reordersam.nf: Memory request adjustment. Adjust PICARD temp directory to Nextflow work directory.  
18. picard/picard_sortsam.nf: Wallclock request adjustment. 
19. python/python_lymphoma_classifier.nf: Typo correction in output name.
20. python/python_somatic_vcf_finalization.nf: Added explicit genome support to facilitate adding mouse to PTA.
21. python/python_split_mnv.nf: Memory request adjustment.
22. r/annotate_sv.nf: Added explicit genome support to facilitate adding mouse to PTA.
23. r/annotate_sv_with_cnv.nf: Minor output file name adjustment.
24. rsem/rsem_alignment_expression.nf: Memory request adjustment. Remove dynamic memory request for STAR genome sort to correct memory failure errors. Added support to save STAR alignment logs.
25. samtools/samtools_filter_unique_reads.nf: Adjust expected file name input.
26. snpeff_snpsift/snpsift_annotate.nf: Adjusted output file name with respect to PTA.
27. svaba/svaba.nf: Adjust Nextflow output streams to caputure index files.
28. utility_modules/jax_trimmer.nf: Wallclock request adjustment.
29. xenome/xenome.nf: Wallclock and memory request adjustment. Adjusted temp directory for `fastq-sort` to Nextflow work directory.  
30. All modules: `${task.memory}` replaced the incorrect `${task.mem}` in the Nextflow error catch statement. 

### Script Added: 

1. pta/annotate-bedpe-with-genes-mouse.r: Removed human specific database expectations. 
2. pta/annotate-cnv-delly.r: Adjusted CNV annotation for Delly output.
3. pta/delly_cnv_plot.r: Added Delly CNV plot. 

### Script Changes: 

1. pta/annotate-bedpe-with-databases.r: Added genome support. For BED annotations, the existing script checks for ANY overlap between BED intervals. For mouse data, this lead to errant overlaps in small InDEL and inversion regions; therefore, mouse PTA requires 80% overlap between target region and query BED.
2. pta/filter-bedpe.r: For mouse PTA we know the type of SV event annotated from databases; therefore, we filter only calls that match annotation type (i.e., DEL, INS, INV). Adjustment to CNV breakpoint checks for cases when breakpoints are not present for targets being annotated. This can occur in mouse PTA due to the change to Delly CNV calling.
4. pta/make_main_vcf.py: Added explicit genome support to facilitate adding mouse to PTA.
5. pta/make_txt.py: Added explicit genome support to facilitate adding mouse to PTA.
6. pta/merge-caller-vcfs.r: Added support for Delly. For Manta the 'infer missing breakpoint' was added as the caller does not insert the reciprocal call in the VCF as the other callers do.

## Release 0.4.5

In this minor release we have updated GBRS and EMASE containers to include a correction made on an index position bug in GBRS genotype printing. GBRS was failing to print the final gene genotype on each chromosome to the `*.genotype.tsv` file.

### Pipelines Added:

None

### Modules Added:

None

### Pipeline Changes:

None

### Module Changes:

1. All EMASE and GBRS modules updated to the latest version of the EMASE/GBRS container. 


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