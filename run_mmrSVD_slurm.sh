#!/bin/bash
#SBATCH --job-name=nf_mmrSVD
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 7-00:00:00
#SBATCH --mem=2000

export NXF_SINGULARITY_CACHEDIR=/path/to/sing_cache
	
module load singularity
	
nextflow -c /path/to/params.config run /path/to/mmrSVD/main.nf -profile slurm,singularity --genome mm10