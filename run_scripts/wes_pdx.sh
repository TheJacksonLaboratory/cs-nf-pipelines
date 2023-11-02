#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=wes_pdx_human
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=1G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow

# RUN PIPELINE
nextflow ../main.nf \
--workflow pdx_wes \
-profile sumner \
--sample_folder <PATH_TO_YOUR_SEQUENCES> \
--gen_org human \
--pubdir "/flashscratch/${USER}/outputDir" \
-w "/flashscratch/${USER}/outputDir/work" \
--comment "This script will run whole exome sequencing analysis on pdx samples using default hg38"
