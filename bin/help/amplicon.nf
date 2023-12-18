def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--organize_by | sample | How to organize the output folder structure. Options: sample or analysis.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.
--sample_folder | /<PATH> | The path to the folder that contains all the samples to be run by the pipeline. The files in this path can also be symbolic links. 
--extension | .fastq.gz | The expected extension for the input read files.
--pattern | '*_R{1,2}*' | The expected R1 / R2 matching pattern. The default value will match reads with names like this READ_NAME_R1_MoreText.fastq.gz or READ_NAME_R1.fastq.gz
--read_type | PE | Options: PE and SE. Default: PE. Type of reads: paired end (PE) or single end (SE).
--concat_lanes | false | Options: false and true. Default: false. If this boolean is specified, FASTQ files will be concatenated by sample. This option is used in cases where samples are divided across individual sequencing lanes.
--csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.
--download_data | null | Requires `--csv_input`. When specified, read data in the CSV manifest will be downloaded from provided URLs. 
--multiqc_config | /<PATH> | The path to amplicon.yaml. The configuration file used while running MultiQC
--cutadaptMinLength | 20 | The minimum length to discard processed reads.
--cutadaptQualCutoff | 20 | The quality cutoff used to trim low-quality ends from reads.
--cutadaptAdapterR1 | 'AGATCGGAAGAG' | TruSeq Illumina adapter trimmer sequence for cutadapt. Change to sequence required by sample libraries. 
--cutadaptAdapterR2 | 'AGATCGGAAGAG' | TruSeq Illumina adapter trimmer sequence. Change to sequence required by sample libraries.  
--ref_fa | '/projects/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta' | The reference fasta to be used throughout the process for alignment as well as any downstream analysis. JAX users should not change this parameter.
--ref_fa_indices | '/projects/omics_share/human/GRCh38/genome/indices/gatk/bwa/Homo_sapiens_assembly38.fasta' | Pre-compiled BWA index files. JAX users should not change this parameter.
--mismatch_penalty | -B 8 | The BWA penalty for a mismatch.
--masterfile | '/projects/compsci/omics_share/human/GRCh38/supporting_files/capture_kit_files/IDT/xGen_sampleID_amplicon/hg38Lifted_xGen_masterfile.txt' | Primerclip master file for amplicon primer trimming. This file is specific to each IDT xGen amplicon panel, and should be changed if samples are not derived from xGen sample ID. 
--amplicon_primer_intervals | '/projects/compsci/omics_share/human/GRCh38/supporting_files/capture_kit_files/IDT/xGen_sampleID_amplicon/hg38Lifted_xGen_SampleID_primers.interval_list' | GATK interval file for calculation of coverage metrics. This file is specific to each IDT xGen amplicon panel, and should be changed if samples are not derived from xGen sample ID.  
--amplicon_target_intervals | '/projects/compsci/omics_share/human/GRCh38/supporting_files/capture_kit_files/IDT/xGen_sampleID_amplicon/hg38Lifted_xGen_SampleID_merged_targets.interval_list' | GATK interval file for calculation of coverage metrics. This file is specific to each IDT xGen amplicon panel, and should be changed if samples are not derived from xGen sample ID.  
--amplicon_rsid_targets | '/projects/compsci/omics_share/human/GRCh38/supporting_files/capture_kit_files/IDT/xGen_sampleID_amplicon/hg38Lifted_xGen_SampleID_merged_targets.txt' | Amplicon SNP target file containing rsID and gene information. Used in generation of the final fingerprint report file.  
--gold_std_indels | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz’ | Human Only - Used in GATK BaseRecalibrator. JAX users should not change this parameter.
--phase1_1000G | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/1000G_phase1.snps.high_confidence.hg38.vcf.gz' | Human Only - Used in GATK BaseRecalibrator. JAX users should not change this parameter.
--dbSNP | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz' | The dbSNP database contains known single nucleotide polymorphisms, and is used in the annotation of known variants.
--dbSNP_index | '/projects/omics_share/human/GRCh38/genome/annotation/snps_indels/dbsnp_151.vcf.gz.tbi' | The dbSNP index file associated with the dbSNP VCF file. 
--call_val | 50 | The minimum phred-scaled confidence threshold at which variants should be called.
--ploidy_val | '-ploidy 2' | variable in haplotypecaller. Not required for amplicon, but present in module. 
--target_gatk | '/projects/compsci/omics_share/human/GRCh38/supporting_files/capture_kit_files/IDT/xGen_sampleID_amplicon/hg38Lifted_xGen_SampleID_merged_targets.bed' | A bed file with amplicon target intervals as defined in the amplicon array used in the data. NOTE: This file MUST reflect the amplicon array used to generate your data.
--call_val | 50 | The minimum phred-scaled confidence threshold at which variants should be called.
'''
}
