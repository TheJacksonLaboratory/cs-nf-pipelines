#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=atac_human
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=1G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/23.10.1

# RUN PIPELINE
nextflow ../main.nf \
--workflow atac \
-profile sumner2 \
--sample_folder <PATH_TO_YOUR_SEQUENCES> \
--gen_org human \
--genome_build 'GRCh38' \
--pubdir "/flashscratch/${USER}/outputDir" \
-w "/flashscratch/${USER}/outputDir/work" \
--comment "This script will run atac sequencing analysis on human samples using default hg38"
