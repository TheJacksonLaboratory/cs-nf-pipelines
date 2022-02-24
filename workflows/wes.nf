#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/wes.nf'
include {param_log} from '../bin/log/wes.nf'
include {BWA_MEM} from '../modules/bwa'
include {SAMTOOLS_INDEX} from '../modules/samtools'
include {READ_GROUPS} from '../modules/read_groups'
include {QUALITY_STATISTICS} from '../modules/quality_stats'
include {AGGREGATE_STATS} from '../bin/wes/aggregate_stats_wes'
include {CAT_HUMAN;
         CAT_HUMAN as CAT_HUMAN_SNP;
         CAT_HUMAN as CAT_HUMAN_INDEL} from '../bin/wes/cat'
include {COSMIC_ANNOTATION;
        COSMIC_ANNOTATION as COSMIC_ANNOTATION_SNP;
        COSMIC_ANNOTATION as COSMIC_ANNOTATION_INDEL} from '../modules/cosmic'
include {PICARD_SORTSAM;
         PICARD_MARKDUPLICATES;
         PICARD_COLLECTHSMETRICS} from '../modules/picard'
include {SNPEFF;
         SNPEFF_HUMAN as SNPEFF_HUMAN_SNP;
         SNPEFF_HUMAN as SNPEFF_HUMAN_INDEL} from '../modules/snpeff'
include {SNPSIFT_EXTRACTFIELDS;
         SNPSIFT_EXTRACTFIELDS as SNPSIFT_EXTRACTFIELDS_SNP;
         SNPSIFT_EXTRACTFIELDS as SNPSIFT_EXTRACTFIELDS_INDEL;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from '../modules/snpsift'
include {GATK_HAPLOTYPECALLER;
         GATK_HAPLOTYPECALLER as GATK_HAPLOTYPECALLER_GVCF;
         GATK_INDEXFEATUREFILE;
         GATK_VARIANTFILTRATION;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL;
         GATK_VARIANTANNOTATOR;
         GATK_MERGEVCF;
         GATK_SELECTVARIANTS;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL;
         GATK_BASERECALIBRATOR;
         GATK_APPLYBQSR} from '../modules/gatk'

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

// main workflow
workflow WES {
  // Step 1: Qual_Stat
  QUALITY_STATISTICS(read_ch)
  // Step 2: Get Read Group Information
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq)
  // Step 3: BWA-MEM Alignment
  BWA_MEM(QUALITY_STATISTICS.out.trimmed_fastq, READ_GROUPS.out.read_groups )
  // Step 4: Variant Preprocessing - Part 1
  PICARD_SORTSAM(BWA_MEM.out.sam)
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)

  // If Human: Step 5-10
  if (params.gen_org=='human'){

    // Step 5: Variant Pre-Processing - Part 2
      GATK_BASERECALIBRATOR(PICARD_MARKDUPLICATES.out.dedup_bam)
      GATK_APPLYBQSR(PICARD_MARKDUPLICATES.out.dedup_bam,
                     GATK_BASERECALIBRATOR.out.table)

    // Step 6: Variant Pre-Processing - Part 3
      PICARD_COLLECTHSMETRICS(GATK_APPLYBQSR.out.bam,
                              GATK_APPLYBQSR.out.bai)

    // Step 7: Variant Calling
      GATK_HAPLOTYPECALLER(GATK_APPLYBQSR.out.bam,
                           GATK_APPLYBQSR.out.bai,
                          'variant')
      GATK_HAPLOTYPECALLER_GVCF(GATK_APPLYBQSR.out.bam,
                                GATK_APPLYBQSR.out.bai,
                               'gvcf')

    // Step 8: Variant Filtration
      // SNP
        GATK_SELECTVARIANTS_SNP(GATK_HAPLOTYPECALLER.out.vcf,
                                GATK_HAPLOTYPECALLER.out.idx,
                               'SNP')
        GATK_VARIANTFILTRATION_SNP(GATK_SELECTVARIANTS_SNP.out.vcf,
                               GATK_SELECTVARIANTS_SNP.out.idx,
                              'SNP')

      // INDEL
      	GATK_SELECTVARIANTS_INDEL(GATK_HAPLOTYPECALLER.out.vcf,
                                  GATK_HAPLOTYPECALLER.out.idx,
                                 'INDEL')
        GATK_VARIANTFILTRATION_INDEL(GATK_SELECTVARIANTS_INDEL.out.vcf,
                                     GATK_SELECTVARIANTS_INDEL.out.idx,
                                    'INDEL')

    // Step 9: Post Variant Calling Processing - Part 1
      // SNP
        COSMIC_ANNOTATION_SNP(GATK_VARIANTFILTRATION_SNP.out.vcf)
        SNPEFF_HUMAN_SNP(COSMIC_ANNOTATION_SNP.out.vcf, 'SNP')
        SNPSIFT_DBNSFP_SNP(SNPEFF_HUMAN_SNP.out.vcf, 'SNP')
        CAT_HUMAN_SNP(SNPSIFT_DBNSFP_SNP.out.vcf, 'SNP')
        SNPSIFT_EXTRACTFIELDS_SNP(CAT_HUMAN_SNP.out.vcf)

      // INDEL
        COSMIC_ANNOTATION_INDEL(GATK_VARIANTFILTRATION_INDEL.out.vcf)
        SNPEFF_HUMAN_INDEL(COSMIC_ANNOTATION_INDEL.out.vcf, 'INDEL')
        SNPSIFT_DBNSFP_INDEL(SNPEFF_HUMAN_INDEL.out.vcf, 'INDEL')
        CAT_HUMAN_INDEL(SNPSIFT_DBNSFP_INDEL.out.vcf, 'INDEL')
        SNPSIFT_EXTRACTFIELDS_INDEL(CAT_HUMAN_INDEL.out.vcf)

    // Step 10: Post Variant Calling Processing - Part 2
      GATK_MERGEVCF(CAT_HUMAN_SNP.out.vcf,
                    CAT_HUMAN_INDEL.out.vcf)


  } else if (params.gen_org=='mouse'){

    // Step 6: Variant Pre-Processing - Part 3
      PICARD_COLLECTHSMETRICS(PICARD_MARKDUPLICATES.out.dedup_bam,
                              PICARD_MARKDUPLICATES.out.dedup_bai)

    // Step 7: Variant Calling
      GATK_HAPLOTYPECALLER(PICARD_MARKDUPLICATES.out.dedup_bam,
                           PICARD_MARKDUPLICATES.out.dedup_bai,
                          'varient')

    // Step 8: Variant Filtration
      GATK_VARIANTFILTRATION(GATK_HAPLOTYPECALLER.out.vcf,
                             GATK_HAPLOTYPECALLER.out.idx,
                            'BOTH')

    // Step 10: Post Variant Calling Processing
      SNPEFF(GATK_VARIANTFILTRATION.out.vcf)
      GATK_VARIANTANNOTATOR(GATK_VARIANTFILTRATION.out.vcf,
                            SNPEFF.out.vcf)
      SNPSIFT_EXTRACTFIELDS(GATK_VARIANTFILTRATION.out.vcf)

  }

  // Step 11: Aggregate Stats
  AGGREGATE_STATS(QUALITY_STATISTICS.out.quality_stats,
						      PICARD_COLLECTHSMETRICS.out.hsmetrics,
						      PICARD_MARKDUPLICATES.out.dedup_metrics)
}
