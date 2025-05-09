#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=atac_mouse
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
--workflow atac \
-profile sumner2 \
--sample_folder <PATH_TO_YOUR_SEQUENCES> \
--gen_org mouse \
--genome_build 'GRCm38' \
--effective_genome_size 2652783500 \
--bowtie2Index '/projects/omics_share/mouse/GRCm38/genome/indices/ensembl/v102/bowtie2/Mus_musculus.GRCm38.dna.primary_assembly.fa' \
--chain '' \
--pubdir "/flashscratch/${USER}/outputDir" \
-w "/flashscratch/${USER}/outputDir/work" \
--comment "This script will run atac sequencing analysis on mouse samples using default mm10"
