#!/bin/bash
#SBATCH --job-name=nf_hg38_RNA_Exp_Est_PE_stranded
#SBATCH --mail-type=END
#SBATCH --mail-user=first.last@jax.org
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 48:00:00
#SBATCH --mem=2000

START_TIME=$(date +'%Y%m%d.%H.%M')     # get date time to be used in default output folder

# original name of file 'run_hg38_RNASeq_PE_stranded.sh'
# HELP FOR USE OF THIS SCRIPT
Help()
{
   # Display Help
   echo
   echo "###############################################################################"
   echo
   echo "RNASeq Pipeline Nextflow and Singularity."
   echo "Authors: Carolyn Paisie, Brian Sanderson, Mike Lloyd, & Barry Guglielmo"
   echo
   echo "###############################################################################"
   echo
   echo "Simplist Usage:"
   echo "   sbatch ./run_rnaseq.sh -f /path/to/samples"
   echo
   echo "###############################################################################"
   echo
   echo "options:"
   echo '-f  | Fastq Path                    | Required | No Default'
   echo '-o  | Out Directory                 | Optional | Default "./output"'
   echo '-c  | Config File                   | Optional | Default "nextflow.config"'
   echo '-t  | Temporary Directory           | Optional | Default "."'
   echo '-m  | Min % High Quality Reads      | Optional | Default "50"'
   echo '-s  | Seed Length RSEM              | Optional | Default "25"'
   echo '-e  | File Extension                | Optional | Default ".fastq.gz"'
   echo '-rp | Type of RNA Prep              | Optional | Default "stranded"'
   echo '         options: "stranded", "..."'
   echo '-rt | Type of Sequencing Reads      | Optional | Default "PE"'
   echo '         options: "PE", "..."'
   echo '-sm | Skip Making Config            | Optional | Default "false"'
   echo '-st | Run Test on Sumner            | Optional | Default "false"'
   echo
   echo "NOTE: Reqires Nextflow Installed in Home Directory and for RNASeq.nf in Same Directory as run_rnaseq.sh"
   echo
}

# TAKE IN ALL PARAMITERS TO USE
while getopts f:o:c:t:m:s:e:rp:rt:st:h flag
do
    case "${flag}" in
        f) FQ_PATH=${OPTARG};;           # fastq path
        o) OUTDIR=${OPTARG};;            # out path
        c) CONFIG_FILE=${OPTARG};;       # name of the config file
        t) TMP_DIR=${OPTARG};;           # directory for temporary files
        m) MIN_PCT_HQ_READS=${OPTARG};;  # minimum % high quality reads
        s) SEED_LENGTH=${OPTARG};;       # Seed length for RSEM alignment step
        e) EXTENSION=${OPTARG};;         # File extension
        rp) READ_PREP=${OPTARG};;        # type of RNA preparation protocol
        rt) READS=${OPTARG};;            # type of sequencing reads; paired end
        st) SUMNER_TEST=${OPTARG};;      # option to do test run on sumner
        h) Help;exit 1;;
       \?) echo "Unknown option -$OPTARG";echo "Use -h for help, options, and defaults"; exit 1;;
    esac
done

# SHARED CACHE FOR CONTAINER IMAGES
CACHE='/projects/omics_share/meta/containers/'
# SET VARIABLES TO DEFAULT IF BLANK
if [ -z "${FQ_PATH}" ]; then Help; echo "No FASTQ PATH"; echo "Use -h for help, options and defaults"; exit 1; fi
if [ -z "${OUTDIR}" ]; then OUTDIR='./output';mkdir output || echo "Re-Writing 'output' Folder"; fi
if [ -z "${CONFIG_FILE}" ]; then CONFIG_FILE="nextflow.config"; fi
if [ -z "${TMP_DIR}" ]; then TMP_DIR=$OUTDIR; fi
# COPY CONFIG FILE TO OUTDIR
cp $CONFIG_FILE $OUTDIR/$CONFIG_FILE
# MAKE CHANGES TO CONFIG IF NECESSARY
sed -i '' 's/hey/bob/' $OUTDIR/$CONFIG_FILE
if [ ! -z "${MIN_PCT_HQ_READS}" ]; then
  sed -i '' "s/EXAMPLE=this/EXAMPLE=${MIN_PCT_HQ_READS}/" $OUTDIR/$CONFIG_FILE;
fi
if [ ! -z "${SEED_LENGTH}" ]; then SEED_LENGTH='25'; fi
if [ ! -z "${EXTENSION}" ]; then EXTENSION='.fastq.gz'; fi
if [ ! -z "${READ_PREP}" ]; then READ_PREP='stranded'; fi
if [ ! -z "${READS}" ]; then READS='PE'; fi
if [ ! -z "${SKIP_MAKE}" ]; then SKIP_MAKE='false'; fi
# depth of coverage and such (pdx specific) need to be removed from RNASeq.nf
# variant calling optional -- need to add this in (v2)
# keeping images sepeate (modularize)
# test quick on toy*** example
  # dsl 2
  # point to external containers (Quay.io and Biocontainers) -- asume open source -- no credentials needed as of now
  # default cache = /projects/omics_share/meta/containers/

# RUN NEXTFLOW PIPELINE USING OUR CONFIGS AND FILE LOCATIONS
# ~/nextflow \
# -c ${pc_name}_hg38_params.config \
# run \
# RNASeq.nf \
# -profile slurm,singularity
