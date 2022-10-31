# RELEASE NOTES

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