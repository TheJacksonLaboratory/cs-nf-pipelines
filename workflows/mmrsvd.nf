#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/mmrsvd'
include {param_log} from '../bin/log/mmrsvd'
include {QUALITY_STATISTICS} from '../modules/quality_stats'
//include {NGMLR_MAP} from '../modules/ngmlr'
include {PICARD_SORTSAM;
         PICARD_SORTSAM as PICARD_SORTSAM_ALIGNED;
       	 PICARD_SORTSAM	as PICARD_SORTSAM_DISCORDANT;
       	 PICARD_SORTSAM	as PICARD_SORTSAM_SPLITS;
         PICARD_MARKDUPLICATES} from '../modules/picard'
//include {SNIFFLES} from '../modules/sniffles'
//include {SVIM_ALIGNMENT} from '../modules/svim'
//include {CUTESV_CSS} from '../modules/cutesv'
//include {PBMM2_ALIGN} from '../modules/cutesv'
//include {PBSV_DISCOVERY; PBSV_CALL} from '../modules/pbsv'
include {READ_GROUPS} from '../modules/read_groups'
include {BWA_MEM} from '../modules/bwa'
include {SAMTOOLS_STATS;
         SAMTOOLS_VIEW} from '../modules/samtools'
include {SAMBLASTER} from '../modules/samblaster'
include {LUMPY_EXTRACT_SPLITS;
         LUMPY_CALL_SV} from '../modules/lumpy'
include {BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_LUMPY;
         BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_BREAKDANCER} from '../modules/bcftools'
//include {MANTA_CALLSV} from '../modules/manta'
//include {DELLY_CALL} from '../modules/delly'

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// prepare reads channel
if(params.seq_method=='illumina' && params.extension != '.bam'){
  if (params.read_type == 'PE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
  }
}

if (params.seq_method=='pacbio' && params.extension != '.bam'){
  read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
}

if (params.extension == '.bam'){
  read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
}


// if channel is empty give error message and exit
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

// think about: failing if too many samples and no override (long queue 14 days (batch 3 days max))
// need cleaning step for successful run

// main workflow
workflow MMRSVD {
  /* PACBIO
  if (params.seq_method == 'pacbio') {

     // no quality or trimming on long reads...long reads notoriously lower quality scores!
     // think about adding quality step
     QUALITY_STATISTICS(read_ch)
     NGMLR_MAP(QUALITY_STATISTICS.out.trimmed_fastq)
     PICARD_SORTSAM(NGMLR_MAP.out.sam)
     // Update to v2
     SNIFFLES(PICARD_SORTSAM.out.bam)
     SVIM_ALIGNMENT(PICARD_SORTSAM.out.bam)
     // HANDLE IN CONFIG
     CUTESV(PICARD_SORTSAM.out.bam,
            PICARD_SORTSAM.out.bai)

     // moved pbmm2 index (.mmi) to be a params instead of needing to generate here for each sample
     // MOVE THE PARAMS TO CONFIG
     PBMM2_ALIGN(QUALITY_STATISTICS.out.trimmed_fastq)

     // ccs and clr command differs and is handled in module
     // CONFIG discovery defult = /ref_data/mm10_ucsc_simple_tandem_repeats_ucsc.bed
     PBSV_DISCOVERY(PBMM2_ALIGN.out.bam)

     // ccs and clr differ only in "--ccs" needed for ccs (Good)
     // tandem vs no tandemm handled in module
     PBSV_CALL(PBSV_DISCOVERY.out.svsig)
     SURVIVOR_MERGE([List of VCFs Here])
  }
*/

  // ILLUMINA
  if (params.seq_method == 'illumina') {
    // if input is raw reads
    if (params.extension != ".bam"){
    	QUALITY_STATISTICS(read_ch)
    	READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq, "")
    	BWA_MEM(QUALITY_STATISTICS.out.trimmed_fastq,
    	        READ_GROUPS.out.read_groups)
    	PICARD_SORTSAM(BWA_MEM.out.sam, "")
        PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
     }
 
    // if input is mapped reads
    // is it okay to assume this is where .bam analysis starts?
    if (params.extension == '.bam') {
  	 PICARD_MARKDUPLICATES(read_ch)
    }

    // mapped or raw continue with same steps from here
    SAMTOOLS_STATS(PICARD_MARKDUPLICATES.out.dedup_bam)
    SAMBLASTER(PICARD_MARKDUPLICATES.out.dedup_bam,
               PICARD_MARKDUPLICATES.out.dedup_bai)
    SAMTOOLS_VIEW(PICARD_MARKDUPLICATES.out.dedup_bam)

    // may need to add naming option for these sort steps
    // is there a reason they would be out of order now? Is there threading going on above?
    PICARD_SORTSAM_ALIGNED(SAMBLASTER.out.bam, "_aligned")
    PICARD_SORTSAM_DISCORDANT(SAMTOOLS_VIEW.out.bam, "_discordant")
    LUMPY_EXTRACT_SPLITS(PICARD_SORTSAM_ALIGNED.out.bam)
    // Keep It (make sure name sorting [etc.], indexing reasons as well)
    PICARD_SORTSAM_SPLITS(LUMPY_EXTRACT_SPLITS.out.bam, "_splits")
    // also includes vcfSort.sh being called. Should we move that call?
    LUMPY_CALL_SV(PICARD_SORTSAM_ALIGNED.out.bam,
                  PICARD_SORTSAM_DISCORDANT.out.bam,
                  PICARD_SORTSAM_SPLITS.out.bam)
    // VCFSORT_SH() // FUTURE

    // COLLECT VCFS THEN REHEAD TO SAMPLEID + CALLER
    // MANTA REHEADER WILL NOT WORK -- MANTA SPECIAL PROCESS
    BCFTOOLS_REHEADER_LUMPY(LUMPY_CALL_SV.out.vcf, 'LUMPY')
/*
    // requires bam2cfg.pl for formatting.
    // Break into own module.
    BAM2CFG(PICARD_MARKDUPLICATES.out.bam)
  	BREAKDANCER_MAX(BAM2CFG.out.bam)

    // breakdancer2vcfHeader.py followed by vcfSort.sh
  	FORMAT_BREAKDANCER(BREAKDANCER_MAX.out.vcf)

  	BCFTOOLS_REHEADER_BREAKDANCER(FORMAT_BREAKDANCER.out.vcf)

    MANTA_CALLSV(PICARD_MARKDUPLICATES.out.bam,
                 PICARD_MARKDUPLICATES.out.bai)
    DELLY_CALL(PICARD_MARKDUPLICATES.out.bam,
               PICARD_MARKDUPLICATES.out.bai)

    // Necessary? BCFTOOLS_VIEW() | vcfSort.sh | reheader

    SURVIVOR_MERGE([List of VCFs Here])
*/  
  }
}

  /* Continue pipeline samd for long or short reads
  // Simple bash script to run surv_annot.sh
  ANNOTATE_SV(SURVIVOR_MERGE.out.vcf)

  // sv_to_table.py
  SUMMARIZE_SV(SURVIVOR_MERGE.out.vcf)

  // RSCRIPT
  PREP_BEDS(ANNOTATE_SV.out.txt,
            SUMMARIZE_SV.out.csv)

  // make the out names a little more descriptive?
  INTERSECT_BEDS(PREP_BEDS.out.ins,
                 PREP_BEDS.out.inv,
                 PREP_BEDS.out.del,
                 PREP_BEDS.out.dup,
                 PREP_BEDS.out.tra)

  // Is this expecting exact names to find? or just ALL .bed files?
  // may be a less explicit way of doing this, but if this function is
  // unique to this pipeline we wont worry about it
  SUMMARIZE_INTERSECTIONS(SUMMARIZE_SV.out.csv,
                          PREP_BEDS.out.ins,
                          PREP_BEDS.out.inv,
                          PREP_BEDS.out.del,
                          PREP_BEDS.out.dup,
                          PREP_BEDS.out.tra,
                          INTERSECT_BEDS.out.ins_s,
                          INTERSECT_BEDS.out.ins_e,
                          INTERSECT_BEDS.out.del_s,
                          INTERSECT_BEDS.out.del_e,
                          INTERSECT_BEDS.out.inv_e,
                          INTERSECT_BEDS.out.tra_e,
                          INTERSECT_BEDS.out.dup_e,
                          INTERSECT_BEDS.out.ins_genes,
                          INTERSECT_BEDS.out.del_genes,
                          INTERSECT_BEDS.out.inv_genes,
                          INTERSECT_BEDS.out.dup_genes,
                          INTERSECT_BEDS.out.tra_genes,
                          INTERSECT_BEDS.out.ins_exons,
                          INTERSECT_BEDS.out.del_exons,
                          INTERSECT_BEDS.out.inv_exons,
                          INTERSECT_BEDS.out.dup_exons,
                          INTERSECT_BEDS.out.tra_exons)

  // Python file
  ANNOTATE_EXONS(INTERSECT_BEDS.out.ins_exons,
                 INTERSECT_BEDS.out.del_exons,
                 INTERSECT_BEDS.out.inv_exons,
                 INTERSECT_BEDS.out.dup_exons,
                 INTERSECT_BEDS.out.tra_exons)

}
// NOOS-33-mmrsvd-pipeline-to-dsl2
*/
