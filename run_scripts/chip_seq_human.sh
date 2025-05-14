#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=rnaseq_human
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=5G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/24.10.6

# RUN PIPELINE
nextflow ../main.nf \
--workflow chipseq \
-profile sumner2 \
--gen_org human \
--genome_build 'GRCh38' \
--input '<PATH_TO_YOUR_SEQUENCES/CSV_input.csv' \
--pubdir "/flashscratch/${USER}/outputDir" \
-w "/flashscratch/${USER}/outputDir/work" \
--narrow_peak \
--comment "This script will run ChIP-sequencing analysis on human samples using default hg38"
