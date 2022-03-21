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

     // no quality or trimming on long reads...long reads notoriously lower quality scores!
     // think about adding quality step
  	 QUALITY_STATISTICS(read_ch)

     // why only fq1 used?
  	 NGMLR_MAP(QUALITY_STATISTICS.out.trimmed_fastq)

     // swapped to picard sortsam rather than samtools sort. is that okay?
     PICARD_SORTSAM(NGMLR_MAP.out.sam)

     // sniffles up to v2 (good to update?)
     SNIFFLES(NGMLR_SORT.out.bam)
     SVIM_ALIGNMENT(NGMLR_SORT.out.bam)

     // ccs and clr differ only in params (this will be handled in config for this step)
  	 CUTESV_CSS(PICARD_SORTSAM.out.bam,
  	 						PICARD_SORTSAM.out.bai)

     // moved pbmm2 index (.mmi) to be a params instead of needing to generate here for each sample
  	 PBMM2_ALIGN(QUALITY_STATISTICS.out.trimmed_fastq,
  	 						 PBMM2_BUILD_INDEX.out.mmi)

     // ccs and clr command differs and is handled in module
  	 PBSV_DISCOVERY(PBMM2_ALIGN.out.bam)

     // ccs and clr differ only in "--ccs" needed for ccs
     // tandem vs no tandemm handled in module
  	 PBSV_CALL(PBSV_DISCOVERY.out.svsig)

     // Add commented code here?
  }

  // ILLUMINA
  if (params.seq_method == 'illumina') {
    // if input is raw reads
    if (params.extension != ".bam"){
    	QUALITY_STATISTICS(read_ch)
    	READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq)
    	BWA_MEM(QUALITY_STATISTICS.out.trimmed_fastq,
    	        READ_GROUPS.out.read_groups)
    	PICARD_SORTSAM(BWA_MEM.out.sam)
      PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
    }

    // if input is mapped reads
    // is it okay to assume this is where .bam analysis starts?
    if (params.extension == '.bam') {
  	 PICARD_MARKDUPLICATES(read_ch)
    }

    // mapped or raw continue with same steps from here
  	// AVOID_RACE_COND() // NULL JUST ECHOS
  	SAMTOOLS_STATS(PICARD_MARKDUPLICATES.out.dedup_bam)
  	SAMBLASTER(PICARD_MARKDUPLICATES.out.dedup_bam
               PICARD_MARKDUPLICATES.out.dedup_bai)
    SAMTOOLS_VIEW_DISCORDANT(PICARD_MARKDUPLICATES.out.dedup_bam
                             PICARD_MARKDUPLICATES.out.dedup_bai)

    // may need to add naming option for these sort steps
    // is there a reason they would be out of order now? Is there threading going on above?
  	PICARD_SORTSAM_ALIGNED(SAMBLASTER.out.bam)
  	PICARD_SORTSAM_DISCORDANT(SAMTOOLS_VIEW_DISCORDANT.out.bam)

    // What exactly is extractSplitReads_BwaMem
  	LUMPY_EXTRACT_SPLITS(LUMPY_MAP.out.aligned)

    // Once again, is it necessary to sort again?
  	PICARD_SORTSAM_SPLITS(LUMPY_EXTRACT_SPLITS.out.bam)

    // the only actual "lumpy" call
    // also includes vcfSort.sh being called. Should we move that call?
  	LUMPY_CALL_SV(PICARD_SORTSAM_ALIGNED.out.bam,
                  PICARD_SORTSAM_DISCORDANT.out.bam,
                  PICARD_SORTSAM_SPLITS.out.bam)

    // are we expecting a single sample per vcf? does breakdancer require spcecial id?
    BCFTOOLS_REHEADER(LUMPY_CALL_SV.out.vcf)

    // requires bam2cfg.pl for formatting.
    // Break into own module.
    BAM2CFG(PICARD_MARKDUPLICATES.out.bam)
  	BREAKDANCER_MAX(BAM2CFG.out.bam)

    // breakdancer2vcfHeader.py followed by vcfSort.sh
  	FORMAT_BREAKDANCER(BREAKDANCER_MAX.out.vcf)

    // is this really neccessary? does the information get lost from above?
    // if we assume one sample per file as default we can just rename at the end.
  	BCFTOOLS_REHEADER(FORMAT_BREAKDANCER.out.vcf)

    MANTA_CALLSV(PICARD_MARKDUPLICATES.out.bam,
                 PICARD_MARKDUPLICATES.out.bai)
    DELLY_CALL(PICARD_MARKDUPLICATES.out.bam,
               PICARD_MARKDUPLICATES.out.bai)

    // Necessary? BCFTOOLS_VIEW() | vcfSort.sh | reheader
  }

  // CONTINUE THE REST OF THE PIPELINE
  SURVIVOR_MERGE()
  ANNOTATE_SV(SURVIVOR_MERGE.out.vcf)
  SUMMARIZE_SV(SURVIVOR_MERGE.out.vcf)
  PREP_BEDS() // RSCRIPT
  INTERSECT_BEDS()
  SUMMARIZE_INTERSECTIONS()
  ANNOTATE_EXONS()

}
