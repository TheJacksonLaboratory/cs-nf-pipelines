#!/bin/bash

#SBATCH --job-name=germline_sv
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=2G
#SBATCH --ntasks=1

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow


export NXF_ANSI_SUMMARY=true
export NXF_ANSI_LOG=true

# RUN TEST PIPELINE
nextflow ../main.nf \
  -profile sumner \
  --workflow germline_sv \
  --data_type illumina \
  --pubdir "/flashscratch/${USER}/mmrsvd_germline" \
  -w "/flashscratch/${USER}/work" \
  --genome_build 'GRCm39' \
  --fastq1 /projects/compsci/omics_share/meta/benchmarking/sim_sv/GRCm39/illumina/GRCm39.v102.autosomes_simSV_germline.R1.fq.gz \
  --fastq2 /projects/compsci/omics_share/meta/benchmarking/sim_sv/GRCm39/illumina/GRCm39.v102.autosomes_simSV_germline.R2.fq.gz  \
  --fasta /projects/omics_share/mouse/GRCm39/genome/sequence/ensembl/v105/Mus_musculus.GRCm39.dna.primary_assembly.fa \
  --fasta_index /projects/omics_share/mouse/GRCm39/genome/sequence/ensembl/v105/Mus_musculus.GRCm39.dna.primary_assembly.fa.fai \
  --bwa_index /projects/omics_share/mouse/GRCm39/genome/indices/ensembl/v105/bwa/Mus_musculus.GRCm39.dna.primary_assembly.fa \
  --sampleID mm10 \
  -resume
