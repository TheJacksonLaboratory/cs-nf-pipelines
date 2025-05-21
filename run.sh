#!/bin/bash

#SBATCH --job-name=CS_nextflow_example
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=2G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/24.10.6

# RUN TEST PIPELINE
nextflow main.nf \
-profile sumner2 \
--workflow rnaseq \
--gen_org mouse \
--sample_folder 'test/rna/mouse' \
--pubdir "/flashscratch/${USER}/outputDir" \
-w "/flashscratch/${USER}/outputDir/work"
