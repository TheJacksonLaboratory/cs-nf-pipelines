#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/mmrsvd'
include {param_log} from '../bin/log/mmrsvd'
include {QUALITY_STATISTICS} from '../modules/quality_stats'
include {NGMLR_MAP} from '../modules/ngmlr'
include {PICARD_SORTSAM;PICARD_MARKDUPLICATES} from '../modules/picard'
include {SNIFFLES} from '../modules/sniffles'
include {SVIM_ALIGNMENT} from '../modules/svim'
include {CUTESV_CSS} from '../modules/cutesv'
include {PBMM2_ALIGN} from '../modules/cutesv'
include {PBSV_DISCOVERY; PBSV_CALL} from '../modules/pbsv'
include {READ_GROUPS} from '../modules/read_groups'
include {BWA_MEM} from '../modules/bwa'
include {SAMTOOLS_STATS} from '../modules/samtools'
include {LUMPY_MAP; LUMPY_CALL_SV; LUMPY_EXTRACT_SPLITS} from '../modules/lumpy'
include {BCFTOOLS_VIEW;BCFTOOLS_REHEADER} from '../modules/bcftools'
include {MANTA_CALLSV} from '../modules/manta'
include {DELLY_CALL} from '../modules/delly'



// NOOS-33-mmrsvd-pipeline-to-dsl2
// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// prepare reads channel
if (params.read_type == 'PE'){
  read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
}
else if (params.read_type == 'SE'){
  read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
}

// if channel is empty give error message and exit
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

// main workflow
workflow MMRSVD {
  // PACBIO
  if (params.seq_method == 'pacbio') {
  	 QUALITY_STATISTICS(read_ch)
  	 NGMLR_MAP(QUALITY_STATISTICS.out.trimmed_fastq)
     PICARD_SORTSAM(NGMLR_MAP.out.sam)
     SNIFFLES(NGMLR_SORT.out.bam)
     SVIM_ALIGNMENT(NGMLR_SORT.out.bam)
  	 CUTESV_CSS(NGMLR_SORT.out.bam,
  	 						NGMLR_SORT.out.bai)
  	 PBMM2_ALIGN(QUALITY_STATISTICS.out.trimmed_fastq,
  	 						 PBMM2_BUILD_INDEX.out.mmi)
  	 PBSV_DISCOVERY(PBMM2_ALIGN.out.bam)
  	 PBSV_CALL()
     // lots was scrapped after this, are we trying to recover it?
  }

  // ILLUMINA
  if (params.seq_method == 'illumina') {
  	QUALITY_STATISTICS(read_ch)
  	READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq)
  	BWA_MEM(QUALITY_STATISTICS.out.trimmed_fastq,
  	        READ_GROUPS.out.read_groups)
  	PICARD_SORTSAM(BWA_MEM.out.sam)
  	PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
  	AVOID_RACE_COND() // NULL JUST ECHOS
  	SAMTOOLS_STATS(PICARD_MARKDUPLICATES.out.dedup_bam)
  	LUMPY_MAP(PICARD_MARKDUPLICATES.out.dedup_bam
              PICARD_MARKDUPLICATES.out.dedup_bai)
  	PICARD_SORTSAM_ALIGNED(LUMPY_MAP.out.aligned)
  	PICARD_SORTSAM_DISCORDANT(LUMPY_MAP.out.discordant)
  	LUMPY_EXTRACT_SPLITS(LUMPY_MAP.out.aligned)
  	PICARD_SORTSAM_SPLITS(LUMPY_EXTRACT_SPLITS.out.bam)
  	LUMPY_CALL_SV(PICARD_SORTSAM_ALIGNED.out.bam,
                  PICARD_SORTSAM_DISCORDANT.out.bam,
                  PICARD_SORTSAM_SPLITS.out.bam)
    BCFTOOLS_REHEADER(LUMPY_CALL_SV.out.vcf)
  	BREAKDANCER_MAX(PICARD_MARKDUPLICATES.out.bam,
                    PICARD_MARKDUPLICATES.out.bai)
  	FORMAT_BREAKDANCER(BREAKDANCER_MAX.out.vcf)
  	BCFTOOLS_REHEADER(FORMAT_BREAKDANCER.out.vcf)
    MANTA_CALLSV(out.bam,
                 out.bai)
    DELLY_CALL(out.bam,
               out.bai)
    BCFTOOLS_VIEW()
  	VCF_SORT()
  }
  /*
  // CONTINUE THE REST OF THE PIPELINE
  SURVIVOR_MERGE()
  ANNOTATE_SV(SURVIVOR_MERGE.out.vcf)
  SUMMARIZE_SV(SURVIVOR_MERGE.out.vcf)
  PREP_BEDS() // RSCRIPT
  INTERSECT_BEDS()
  SUMMARIZE_INTERSECTIONS()
  ANNOTATE_EXONS()
  */
}
