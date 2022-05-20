# RELEASE NOTES

## Release 0.2.0

### Pipelines Added:

1. Bisulfite sequencing RRBS - Mouse & Human

### Modules Added: 

NONE

### Pipeline Changes:

NONE

### Module Changes:

NONE

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