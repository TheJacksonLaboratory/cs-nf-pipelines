def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--organize_by | sample | How to organize the output folder structure. Options: sample or analysis.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /fastscratch or other directory with ample storage.

--input | /<PATH> | The path to the design file that contains all the samples to be run by the pipeline. - For design file format, please see : 'https://nf-co.re/chipseq/1.2.1/usage'   
--extension | .fastq.gz | The expected extension for the input read files.
--pattern | '*_R{1,2}*' | The expected R1 / R2 matching pattern. The default value will match reads with names like this READ_NAME_R1_MoreText.fastq.gz or READ_NAME_R1.fastq.gz
--read_type | PE | Options: PE and SE. Default: PE. Type of reads: paired end (PE) or single end (SE).

--gen_org | mouse | Options: mouse and human.
                        
--fragment_size | 200 |  Number of base pairs to extend single-end reads when creating bigWig files (Default: 200)
--fingerprint_bins | 500000 | Number of genomic bins to use when generating the deepTools fingerprint plot. Larger numbers will give a smoother profile, but take longer to run (Default: 500000)                              
--gtf      | The full path to GTF file for annotating peaks and the GTF file should resemble the Ensembl format
--gene_bed | The full path to BED file for genome-wide gene intervals                                                        

--ref_fa         | Mouse: '/projects/omics_share/mouse/GRCm38/genome/sequence/ensembl/v102/Mus_musculus.GRCm38.dna.primary_assembly.fa' 
                 | Human: '/projects/omics_share/human/GRCh38/genome/indices/gatk/bwa/Homo_sapiens_assembly38.fasta'
--ref_fa_indices | Mouse: '/projects/omics_share/mouse/GRCm38/genome/indices/ensembl/v102/bwa/Mus_musculus.GRCm38.dna.primary_assembly.fa' | The default value for mm10. 
                 | Human: '/projects/omics_share/human/GRCh38/genome/indices/gatk/bwa/Homo_sapiens_assembly38.fasta'
                 | Pre-compiled BWA index files, points to human reference when --gen_org human.

--macs_gsize | Effective genome size parameter required by MACS2 | if this parameter is not specified then the MACS2 peak-calling and differential analysis will be skipped                            
--blacklist | The BED format file, if provided, alignments that overlap with the regions in this file will be filtered out | Please see : 'https://sites.google.com/site/anshulkundaje/projects/blacklists'

--non_directional | true | Options: true and false. Selecting this option for non-directional RRBS libraries will screen quality-trimmed sequences for CAA or CGA at the start of the read and, if found, removes the first two base pairs.

--trimLength | 30 | Discard reads that became shorter than length 'INT' because of either quality or adapter trimming. A value of 0 effectively disables this behaviour.
--qualThreshold | 30 | Trim low-quality ends from reads in addition to adapter removal. For RRBS samples, quality trimming will be performed first, and adapter trimming is carried in a second round. Other files are quality and adapter trimmed in a single pass. The algorithm is the same as the one used by BWA (Subtract INT from all qualities; compute partial sums from all indices to the end of the sequence; cut sequence at the index at which the sum is minimal).
--adapOverlap | 1 | Stringency for overlap with adapter sequence required to trim a sequence. Defaults to a very stringent setting of 1, i.e. even a single base pair of overlapping sequence will be trimmed of the 3' end of any read.
--adaptorSeq | 'AGATCGGAAGAGC' | Adapter sequence to be trimmed. This sequence is the standard Illumina adapter sequence.

--mismatch_penalty | '' | The BWA penalty for a mismatch.
--bwa_min_score | false | Don’t output BWA MEM alignments with score lower than this parameter (Default: false)
--keep_dups | false | Duplicate reads are not filtered from alignments (Default: false)
--keep_multi_map | false | Reads mapping to multiple locations in the genome are not filtered from alignments (Default: false)
--bamtools_filter_pe_config | /<PATH> | The path to bamtools_filter_pe.json for paired end (PE)
--bamtools_filter_se_config | /<PATH> | The path to bamtools_filter_se.json for single end (SE)
                            | The configuration file used while running bamtools filter 

--narrow_peak | false | MACS2 is run by default with the --broad flag. Specify this flag to call peaks in narrowPeak mode (Default: false)
--broad_cutoff | 0.1 | Specifies broad cut-off value for MACS2. Only used when --narrow_peak isnt specified (Default: 0.1)

--skip_preseq | false | Skip Preseq
--skip_peak_qc | false | Skip MACS2 peak QC plot generation (Default: false)
--skip_peak_annotation | false | Skip MACS2 peak QC plot generation (Default: false)
--skip_consensus_peaks | false | Skip consensus peak generation, annotation and counting (Default: false)
--skip_diff_analysis | false | Skip differential binding analysis with DESeq2 (Default: false)
--deseq2_vst | false | Use vst transformation instead of rlog with DESeq2. (Default: false)
--macs_fdr | false | Minimum FDR (q-value) cutoff for peak detection, --macs_fdr and --macs_pvalue are mutually exclusive (Default: false)
--macs_pvalue | false | p-value cutoff for peak detection (Default: false).
--min_reps_consensus | 1 | Number of biological replicates required from a given condition for a peak to contribute to a consensus peak (Default: 1)
--save_macs_pileup | false | Instruct MACS2 to create bedGraph files using the -B --SPMR parameters (Default: false).

--tmpdir  | /<PATH> | Temporary directory to store temp files.  
'''
}

