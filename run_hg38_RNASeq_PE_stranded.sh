#!/bin/bash
#SBATCH --job-name=nf_hg38_RNA_Exp_Est_PE_stranded
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@jax.org
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 48:00:00
#SBATCH --mem=2000



~/nextflow \
-c ${pc_name}_hg38_params.config \
run \
/projects/compsci/nextflow/pipelines/RNA/RNA_Expression_Estimation_Single_Sample/1.0.0/RNASeq.nf \
-profile slurm,singularity



