#!/bin/bash

#SBATCH --job-name=CS_nextflow_example
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 20:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD SINGULARITY
module load singularity

# RUN TEST PIPELINE
~/nextflow main.nf --workflow rna --sample_folder 'test/rna' --pubdir '/fastscratch/outputDir' -w '/fastscratch/outputDir/work'