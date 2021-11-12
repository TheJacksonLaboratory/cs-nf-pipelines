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
   echo "Authors: Brian Sanderson, Mike Lloyd, & Barry Guglielmo"
   echo
   echo "###############################################################################"
   echo
   echo "options:"
   echo '-f  | Fastq Path                    | Required | No Default'
   echo '-o  | Out Directory                 | Optional | Default "./output"'
   echo '-c  | Config File                   | Optional | Default "hg38_params.config"'
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
   echo "NOTE: Reqires Nextflow Installed in Home Directory."
   echo
}

# TAKE IN ALL PARAMITERS TO USE
while getopts f:o:c:t:m:s:e:rp:rt:sm:st:h flag
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
        sm) SKIP_MAKE=${OPTARG};;        # option to skip making config file
        st) SUMNER_TEST=${OPTARG};;      # option to do test run on sumner
        h) Help;exit 1;;
       \?) echo "Unknown option -$OPTARG";echo "Use -h for help, options, and defaults"; exit 1;;
    esac
done


# sSET VARIABLES TO DEFAULT IF BLANK
if [ -z "${FQ_PATH}" ]; then echo "No FASTQ PATH"; echo "Use -h for help, options and defaults"; exit 1; fi
if [ -z "${OUTDIR}" ]; then OUTDIR='./output';mkdir output || echo "Re-Writing 'output' Folder"; fi
if [ -z "${CONFIG_FILE}" ]; then CONFIG_FILE="hg38_params.config"; fi
if [ -z "${TMP_DIR}" ]; then TMP_DIR=$OUTDIR; fi
if [ -z "${MIN_PCT_HQ_READS}" ]; then MIN_PCT_HQ_READS='50'; fi
if [ -z "${SEED_LENGTH}" ]; then SEED_LENGTH='25'; fi
if [ -z "${EXTENSION}" ]; then EXTENSION='.fastq.gz'; fi
if [ -z "${READ_PREP}" ]; then READ_PREP='stranded'; fi
if [ -z "${READS}" ]; then READS='PE'; fi
if [ -z "${SKIP_MAKE}" ]; then SKIP_MAKE='false'; fi

# MAKE CONFIG FILE(S)
FILE_LIST_R1=`ls ${FQ_PATH}/*_R1*${Extension} > $OUTDIR/r1.txt`
FILE_LIST_R2=`ls ${FQ_PATH}/*_R2*${Extension} > $OUTDIR/r2.txt`
n=`cat $OUTDIR/r1.txt | wc -l`
for ((i=1; i<=n; i++)) # * IS IT REALLY NECESSARY TO HAVE A CONFIG FOR EACH SAMPLE PAIR???
  do
    R1=$(head -n $i $OUTDIR/r1.txt | tail -n 1)
    R2=$(head -n $i $OUTDIR/r2.txt | tail -n 1)
    base=$(basename $R1 .fastq.gz)
    echo " // $CONFIG_FILE : `date`
    params
    {
       _fqPath=\""${FQ_PATH}"\"
       fastqInputs = \"${R1}"",""${R2}"\""
       outdir = \"${OUTDIR}\"
       tmpdir = \"${TMP_DIR}\"
       min_pct_hq_reads = \"${MIN_PCT_HQ_READS}\"
       reads = \"${READS}\"
       read_prep = \"${READ_PREP}\"
       seed_length = \"${SEED_LENGTH}\"
       gen_org = \"human\"
       aligner = \"--bowtie2\"
       rsem_ref_prefix = '/projects/compsci/refdata/Human/hg38/Index_Files/Bowtie2/Homo_sapiens.GRCh38.dna.toplevel_chr_mod_1_22_MT_X_Y'
       ref_fa = '/projects/compsci/refdata/Human/hg38/Index_Files/Bowtie2/Homo_sapiens.GRCh38.dna.toplevel_chr_mod_1_22_MT_X_Y.fa'
       ref_flat = '/projects/compsci/refdata/Human/hg38/Index_Files/Bowtie2/refFlat.txt'
       ribo_intervals = '/projects/compsci/refdata/Human/hg38/Index_Files/Bowtie2/interval_rRNA'
       probes = '/projects/compsci/refdata/Human/agilent/hg38_agilent_SureSelect_V4_pChrM_probes_genename.bed'
       ctp_genes = '/projects/compsci/refdata/Human/agilent/359genes_b38_noheader_withNames.bed'
       cov_calc = '/projects/compsci/nextflow/bin/coveragecalculator.py'
    }"                               > $OUTDIR/${base}_${CONFIG_FILE}

  done
  # for ((i=1; i<=n; i++))
  #
  #   do
  #
  #     R1=$(head -n $i r1.txt | tail -n 1)
  #     base=$(basename $R1 .fastq.gz)
  #     sbatch --export=ALL,pc_name=$base ./run_rnaseq.sh
  #
  #   done
# RUN NEXTFLOW PIPELINE USING OUR CONFIGS AND FILE LOCATIONS
# ~/nextflow \
# -c ${pc_name}_hg38_params.config \
# run \
# RNASeq.nf \
# -profile slurm,singularity
