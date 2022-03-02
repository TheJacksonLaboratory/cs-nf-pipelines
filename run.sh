#!/bin/bash
#SBATCH --mail-user=barry.guglielmo@jax.org
#SBATCH --job-name=min_nextflow
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 20:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=16

cd $SLURM_SUBMIT_DIR

# LOAD SINGULARITY
ml singularity
# RUN TEST PIPELINE
~/nextflow main.nf --workflow rnaseq --gen_org mouse --comment "fixing rnaseq mouse" --ribo_intervals=null
