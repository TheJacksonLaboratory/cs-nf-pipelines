#!/bin/bash

#SBATCH --job-name=CS_nextflow_example
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 01:00:00
#SBATCH --mem=2G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow

# RUN TEST PIPELINE
nextflow main.nf \
-profile sumner \
--workflow rnaseq \
--gen_org mouse \
--sample_folder 'test/rna/mouse' \
--pubdir '/fastscratch/outputDir' \
-w '/fastscratch/outputDir/work'