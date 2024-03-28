#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=pta_human
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
--workflow pta \
-profile sumner2 \
--gen_org mouse \
--csv_input ../test/csv_samplesheets/mm_test_input.csv \
--pubdir "/flashscratch/${USER}/outputDir" \
-w "/flashscratch/${USER}/outputDir/work" \
--comment "This script will run paired tumor analysis on test data"
