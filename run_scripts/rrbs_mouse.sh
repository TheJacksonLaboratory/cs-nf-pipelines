#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=rrbs_human
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
--workflow rrbs \
-profile sumner \
--sample_folder <PATH_TO_YOUR_SEQUENCES> \
--gen_org mouse \
--pubdir '/fastscratch/outputDir' \
-w '/fastscratch/outputDir/work' \
--comment "This script will run the reduced-representation bisulfite sequencing analysis pipeline on mouse samples using default mm10"