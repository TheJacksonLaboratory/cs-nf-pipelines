#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=gbrs_mouse
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
-profile sumner \
--workflow gbrs \
--pubdir "/fastscratch/${USER}/outputDir" \
-w /fastscratch/${USER}/outputDir/work \
--sample_folder <PATH_TO_YOUR_SEQUENCES> \
--sample_generation 42 \
--sample_sex F \
--comment "This script will run gbrs analysis on mouse samples"
