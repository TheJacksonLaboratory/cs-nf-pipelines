#!/bin/bash
#SBATCH --job-name=min_nextflow
#SBATCH --mail-type=END
#SBATCH --mail-user=your.email@jax.org
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 02:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=16

cd $SLURM_SUBMIT_DIR

# LOAD SINGULARITY
ml singularity
~/nextflow ../main.nf
