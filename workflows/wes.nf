#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {BWA_MEM} from '../modules/bwa'
include {QUALITY_STATISTICS} from '../modules/quality_stats'
include {READ_GROUPS} from '../modules/read_groups'
include {PICARD_SORTSAM;PICARD_MARKDUPLICATES} from '../modules/picard'


// prepare reads channel
if (params.read_type == 'PE'){
  read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}_*${params.extension}",checkExists:true )
}
else if (params.read_type == 'SE'){
  read_ch = Channel.fromFilePairs("${params.fq_path}/*${params.extension}",checkExists:true, size:1 )
}

// main workflow
workflow WES {
  // Step 1: Qual_Stat
  QUALITY_STATISTICS(read_ch) 
  // Step 2: Get Read Group Information
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq)
  // Step 3: BWA-MEM Alignment
  BWA_MEM(QUALITY_STATISTICS.out.trimmed_fastq, READ_GROUPS.out.read_groups)
  // Step 4: Variant Preprocessing - Part 1
  PICARD_SORTSAM(BWA_MEM.out.bwa_mem)
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.picard_sortsam_bam)
 
  /* Step 5: Variant Pre-Processing - Part 2
  if (params.gen_org=='human'){
    GATK_BASERECALIBRATOR(PICARD_MARKDUPLICATES.out.bam_dedup)
    GATK_APPLYBQSR(GATK_BASERECALIBRATOR.out.recal_data_table)
    SAMTOOLS_INDEX(GATK_APPLYBQSR.out.bam)
  }
  else if (params.gen_org=='mouse'){
    SAMTOOLS_INDEX(PICARD_MARKDUPLICATES.out.bam_dedup)
  }
  // Step 6: Variant Pre-Processing - Part 3
  PICARD_CALCULATEHSMETRICS(SAMTOOLS_INDEX.out.samtools_index)
  // Step 7: Variant Calling
  GATK_HAPLOTYPECALLER(SAMTOOLS_INDEX.out.samtools_index, 'normal')
  GATK_HAPLOTYPECALLER(SAMTOOLS_INDEX.out.samtools_index, 'gvcf')
  // Step 8-11 Human
  if (params.gen_org=='human'){
    // Step 8: Variant Filtration
      // SNP
      GATK_SELECTVARIANTS()
      GATK_INDEXFEATUREFILE()
      GATK_VARIANTFILTRATION()
      // INDEL
      GATK_SELECTVARIANTS()
      GATK_INDEXFEATUREFILE()
      GATK_VARIANTFILTRATION()
    // Step 9: Post Variant Calling Processing - Part 1 (MAY NEED TO SET SOME VARIABLES TO HOLD PROCESS.OUT)
      // SNP
        COSMIC_ANNOTATION()       // params.cosmic_annot
        SNPEFF()
        SNPSIFT()
        CAT_HUMAN(FILES, 'SNP')
        EXTRACTFIELDS()
      // INDEL
        COSMIC_ANNOTATION()       // params.cosmic_annot
        SNPEFF()
        SNPSIFT_DBNSFP()
        CAT_HUMAN(FILES, 'INDEL')
        SNPSIFT_EXTRACTFIELDS()
    // Step 10: Post Variant Calling Processing - Part 2
      GENOMEANALYSISTK()        //GenomeAnalysisTK
    // Step 11: Aggregate Stats (UNIQUE TO BIN/WES)
      AGGREGATE_STATS()
  }

  // Step 8-11 Mouse
  if (params.gen_org=='mouse'){
    // Step 8: Variant Filtration
      GATK_INDEXFEATUREFILE()
      GATK_VARIANTFILTRATION()
    // Step 9: Post Variant Calling Processing - Part 1
      CAT_SNP_INDEL()
    // Step 10: Post Variant Calling Processing - Part 2
      SNPEFF()                  // snpEff
      GENOMEANALYSISTK()
      EXTRACTFIELDS()
    // Step 11: Aggregate Stats (UNIQUE TO BIN/WES)
      AGGREGATE_STATS()
  }
  */
}
