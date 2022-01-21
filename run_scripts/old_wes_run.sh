#!/bin/bash
#SBATCH --job-name=nf_mm10_WES
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@jax.org
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 48:00:00
#SBATCH --mem=2000




~/nextflow \
-c ${pc_name}_mm10_params.config \
run \
/projects/compsci/nextflow/pipelines/Exome/1.0.0/WholeExome.nf \
-profile slurm,singularity
