#!/bin/bash
#SBATCH --job-name=min_nextflow
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=barry.guglielmo@jax.org
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=16

cd $SLURM_SUBMIT_DIR

# LOAD SINGULARITY
ml singularity

# RUN TEST PIPELINE
~/nextflow main.nf \
--workflow rnaseq \
--fq_path *PATH_TO_YOUR_SEQUENCES* \
--gen_org mouse
