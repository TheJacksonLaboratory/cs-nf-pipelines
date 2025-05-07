# Run Scripts


This directory contains example run scripts for each pipeline by species combination. 

Example: 

```
#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=rnaseq_human
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
--workflow rnaseq \
--sample_folder <PATH_TO_YOUR_SEQUENCES> \
--gen_org human \
--pubdir '/flashscratch/outputDir' \
-w '/flashscratch/outputDir/work' \
--comment "This script will run rnaseq on human samples using default hg38"
```

There are several things a user must change before running these scripts: 

1. `--mail-user=first.last@jax.org` should be set to your email.

2. `<PATH_TO_YOUR_SEQUENCES>` must be changed to point at the location of your data files. NOTE: The script assumes '_R{1,2}*.fastq.gz' is the default matching string for FASTQ files. 
    See Wiki documentation for parameters to adjust if this is not the case. 

3. `--pubdir '/flashscratch/outputDir'` and `-w '/flashscratch/outputDir/work'` should be changed to relevant directories. The `-w` work directory can be quite large. Use of `/flashscratch/` is recommended. 


**NOTE:**  

1. These scripts assume they are being run from within `cs-nf-pipelines/run_scripts`. If they are moved to other locations, specify the absolute path to `main.nf` (e.g., `/home/USERNAME/cs-nf-pipelines/main.nf`)

2. Sample data for each workflow and species are provided in cs-nf-pipelines/test/<DATA-TYPE>