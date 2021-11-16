#!/bin/bash
#SBATCH --job-name=nf_hg38_RNA_Exp_Est_PE_stranded
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@jax.org
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 48:00:00
#SBATCH --mem=2000

# LOAD SINGULARITY
ml singularity

# HELP FOR USE OF THIS SCRIPT
Help()
{
   # Display Help
   echo
   echo "###############################################################################"
   echo
   echo "Description of Script: Bare Minimum to Run Pipeline on Sumner"
   echo "Authors: Barry Guglielmo"
   echo
   echo "###############################################################################"
   echo
   echo "Simplist Usage:"
   echo "   sbatch ./run.sh -f /path/to/samples"
   echo
   echo "###############################################################################"
   echo
   echo "options:"
   echo '-f  | Fastq Path                    | Required | No Default'
   echo '-o  | Out Directory                 | Optional | Default "./output"'
   echo '-c  | Config File                   | Optional | Default "nextflow.config"'
   echo
   echo "NOTE: Reqires Nextflow Installed in Home Directory and for minimum.nf in Same Directory as run.sh"
   echo
}

# TAKE IN ALL PARAMITERS TO USE
while getopts f:o:c:h flag
do
    case "${flag}" in
        f) FQ_PATH=${OPTARG};;           # fastq path
        o) OUTDIR=${OPTARG};;            # out path
        c) CONFIG_FILE=${OPTARG};;       # name of the config file
        h) Help;exit 1;;
       \?) echo "Unknown option -$OPTARG";echo "Use -h for help, options, and defaults"; exit 1;;
    esac
done

# SET INPUT DEFAULTS (THIS COULD BE DONE IN NF SCRIPT, BUT IT READS BETTER HERE)
if [ -z "${FQ_PATH}" ]; then Help; echo "No FASTQ PATH"; echo "Use -h for help, options and defaults"; exit 1; fi
if [ -z "${OUTDIR}" ]; then OUTDIR='./output';mkdir output || echo "Re-Writing 'output' Folder"; fi
if [ -z "${CONFIG_FILE}" ]; then CONFIG_FILE="nextflow.config"; fi

# MAKE COPY NF CONFIG FILE FOR RECORDS
cp $CONFIG_FILE $OUTDIR/$CONFIG_FILE
# RUN NEXTFLOW PIPELINE USING OUR CONFIGS AND FILE LOCATIONS
~/nextflow \
-c $CONFIG_FILE \
run \
minimum.nf \
--outdir $OUTDIR \
--fq_path $FQ_PATH
# -profile slurm,singularity
