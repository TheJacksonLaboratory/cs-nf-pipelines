#!/bin/bash

# Sample command to run script to make config files and launch pipeline:
# sbatch make_hg38_config_run_RNASeq_PE_stranded.sh /projects/compsci/USERS/paisic/hsa_fastq_RNA_seq/PE /projects/compsci/USERS/paisic/hsa_hg38_RNASeq_testing_PE/ /projects/compsci/paisic/hsa_hg38_RNASeq_testing_PE 50 PE stranded 25


#IMPORTANT: snpEff.config needs to be in directory one level above data directory containing snpEff database files!!!!


CONFIG_FILE=hg38_params.config  # name of the config file
FQ_PATH="."                     # path to the fastq directory
OUTDIR="."                      # output directory
TMPDIR="."                      # directory for temporary files
MIN_PCT_HQ_READS='50'           # minimum % high quality reads
READS='PE'                      # type of sequencing reads; paired end
READ_PREP='stranded'            # type of RNA preparation protocol
SEED_LENGTH='25'                # Seed length for RSEM alignment step

if [ $# -eq 0 ]; then
    echo "Useage: $0 fastq_dir [ outdir [tmpdir [min_pct_hq_reads [reads [read_prep [seed_length ]]]]]]"
    echo "defaults: $0 . . 50 PE stranded 25"
    exit 1
fi

FQ_PATH=$1
if [ $# -gt 1 ]; then
    OUTDIR=$2
fi
if [ $# -gt 2 ]; then
    TMPDIR=$3
fi
if [ $# -gt 3 ]; then
    MIN_PCT_HQ_READS=$4
fi
if [ $# -gt 4 ]; then
    READS=$5
fi
if [ $# -gt 5 ]; then
    READ_PREP=$6
fi
if [ $# -gt 6 ]; then
    SEED_LENGTH=$7
fi


FILE_LIST_R1=`ls ${FQ_PATH}/*_R1*.fastq.gz > r1.txt`
FILE_LIST_R2=`ls ${FQ_PATH}/*_R2*.fastq.gz > r2.txt`

n=`cat r1.txt | wc -l`

for ((i=1; i<=n; i++))
  do
    R1=$(head -n $i r1.txt | tail -n 1)
    R2=$(head -n $i r2.txt | tail -n 1)
    base=$(basename $R1 .fastq.gz)

    echo " // $CONFIG_FILE : `date`  	                                                                                                           
    params                                                                                                                                 
    {                                           				                                               
       _fqPath=\""${FQ_PATH}"\"
       fastqInputs = \"${R1}"",""${R2}"\""
       outdir = \"${OUTDIR}\"                                                                                    
       tmpdir = \"${TMPDIR}\"                                                                             
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
    }"                               > ${base}_${CONFIG_FILE}

  done

for ((i=1; i<=n; i++))

  do

    R1=$(head -n $i r1.txt | tail -n 1)
    base=$(basename $R1 .fastq.gz)
    sbatch --export=ALL,pc_name=$base ./run_hg38_RNASeq_PE_stranded.sh

  done
