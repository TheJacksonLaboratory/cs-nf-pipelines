#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=gbrs_mouse
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
--workflow prepare_emase \
--pubdir "/flashscratch/${USER}/outputDir" \
-w /flashscratch/${USER}/outputDir/work \
--genome_file_list "/path/to/genome/A.fa,/path/to/genome/B.fa,..." \
--gtf_file_list "/path/to/gtf/A.gtf,/path/to/gtf/B.gtf,..." \
--haplotype_list "A,B,..." \
--comment "This script will run prepare_emase to generate multiway references based on default parameters"
