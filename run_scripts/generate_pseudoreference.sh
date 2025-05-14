#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=generate_pseudoreference_mouse
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
-profile sumner2 \
--workflow generate_pseudoreference \
--pubdir "/flashscratch/${USER}/outputDir" \
-w /flashscratch/${USER}/outputDir/work \
--sample_folder <PATH_TO_YOUR_SEQUENCES> \
--comment "This script will generate_pseudoreference using default parameters"
