#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=wgs_human
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=1G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD SINGULARITY
module load singularity

# RUN TEST PIPELINE
~/nextflow ../main.nf \
--workflow wgs \
--sample_folder <PATH_TO_YOUR_SEQUENCES> \
--gen_org human \
--pubdir '/fastscratch/outputDir' \
-w '/fastscratch/outputDir/work' \
--comment "This script will run whole genome sequencing on human samples using default hg38"
