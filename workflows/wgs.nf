#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/wes.nf'
include {param_log} from '../bin/log/wes.nf'
include {BWA_MEM} from '../modules/bwa'
include {COSMIC_ANNOTATION as COSMIC_ANNOTATION_SNP;
         COSMIC_ANNOTATION as COSMIC_ANNOTATION_INDEL} from '../modules/cosmic'
include {SNPEFF_HUMAN as SNPEFF_HUMAN_SNP;
         SNPEFF_HUMAN as SNPEFF_HUMAN_INDEL} from '../modules/snpeff'
include {CAT_ANNOTATE as CAT_ANNOTATE_SNP;
         CAT_ANNOTATE as CAT_ANNOTATE_INDEL;
         CAT_ONEPERLINE as CAT_ONEPERLINE_SNP;
         CAT_ONEPERLINE as CAT_ONEPERLINE_INDEL} from '../bin/wgs/cat'
include {SNPEFF} from '../modules/snpeff'
include {SNPSIFT_EXTRACTFIELDS;
         SNPSIFT_EXTRACTFIELDS as SNPSIFT_EXTRACTFIELDS_SNP;
         SNPSIFT_EXTRACTFIELDS as SNPSIFT_EXTRACTFIELDS_INDEL;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from '../modules/snpsift'
include {AGGREGATE_STATS_MOUSE;
         AGGREGATE_STATS_HUMAN} from '../bin/shared/aggregate_stats'
include {READ_GROUPS} from '../modules/read_groups'
include {QUALITY_STATISTICS} from '../modules/quality_stats'
include {PICARD_SORTSAM;
         PICARD_MARKDUPLICATES;
         PICARD_COLLECTALIGNMENTSUMARYMETRICS} from '../modules/picard'
include {GATK_REALIGNERTARGETCREATOR;
         GATK_BASERECALIBRATOR;
         GATK_PRINTREADS;
         GATK_INDELREALIGNER;
         GATK_MERGEVCF;
         GATK_VARIANTANNOTATOR;
         GATK_HAPLOTYPECALLER_WGS;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from '../modules/gatk'

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
workflow WGS {
  // Step 1: Qual_Stat
  QUALITY_STATISTICS(read_ch)
  // Step 2: Get Read Group Information
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq)
  // Step 3: BWA-MEM Alignment
  BWA_MEM(QUALITY_STATISTICS.out.trimmed_fastq, READ_GROUPS.out.read_groups)
  // Step 4: Variant Preprocessing - Part 1
  PICARD_SORTSAM(BWA_MEM.out.sam)
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
  // Step 5
  GATK_REALIGNERTARGETCREATOR(PICARD_MARKDUPLICATES.out.dedup_bam)
  GATK_INDELREALIGNER(PICARD_MARKDUPLICATES.out.dedup_bam,
                        GATK_REALIGNERTARGETCREATOR.out.intervals)
  // If Human
  if (params.gen_org=='human'){

  //    Need help sorting out baserecalibrator and printreads
//    GATK_BASERECALIBRATOR(GATK_INDELREALIGNER.out.bam)
//    ran into issue here -BQSR parameter no longer valid, cannot use old printreads gatk3 and baserecal gatk4 together either
//    GATK_PRINTREADS(GATK_INDELREALIGNER.out.bam,
//                    GATK_BASERECALIBRATOR.out.table)
  //    for now skipping printreads and using indelrealigner output instead

    PICARD_COLLECTALIGNMENTSUMARYMETRICS(GATK_INDELREALIGNER.out.bam)
    GATK_HAPLOTYPECALLER_WGS(GATK_INDELREALIGNER.out.bam,
                             GATK_INDELREALIGNER.out.bai)
  }
  
  // If Mouse
  if (params.gen_org=='mouse'){
    PICARD_COLLECTALIGNMENTSUMARYMETRICS(GATK_INDELREALIGNER.out.bam)
    GATK_HAPLOTYPECALLER_WGS(GATK_INDELREALIGNER.out.bam,
                             GATK_INDELREALIGNER.out.bai)
  }

  // SNP
    GATK_SELECTVARIANTS_SNP(GATK_HAPLOTYPECALLER_WGS.out.vcf,
                            GATK_HAPLOTYPECALLER_WGS.out.idx,
                           'SNP')
    GATK_VARIANTFILTRATION_SNP(GATK_SELECTVARIANTS_SNP.out.vcf,
                               GATK_SELECTVARIANTS_SNP.out.idx,
                              'SNP')
  // INDEL
    GATK_SELECTVARIANTS_INDEL(GATK_HAPLOTYPECALLER_WGS.out.vcf,
                              GATK_HAPLOTYPECALLER_WGS.out.idx,
                             'INDEL')
    GATK_VARIANTFILTRATION_INDEL(GATK_SELECTVARIANTS_INDEL.out.vcf,
                                 GATK_SELECTVARIANTS_INDEL.out.idx,
                                'INDEL')

  // Cat Output to vcf-annotate* 
  // Note: Mouse=[CHROM,FROM,TO,ID] Human=[CHROM,POS,ID,REF,ALT]
    CAT_ANNOTATE_SNP(GATK_VARIANTFILTRATION_SNP.out.vcf)
    CAT_ANNOTATE_INDEL(GATK_VARIANTFILTRATION_INDEL.out.vcf)

// Final Post-Processing Steps Slightly Different for Mouse and Human
  
  // If Human
  if (params.gen_org=='human'){
   
    // SNP
      COSMIC_ANNOTATION_SNP(CAT_ANNOTATE_SNP.out.vcf)
      SNPEFF_HUMAN_SNP(COSMIC_ANNOTATION_SNP.out.vcf, 'SNP')
      SNPSIFT_DBNSFP_SNP(SNPEFF_HUMAN_SNP.out.vcf, 'SNP')
      CAT_ONEPERLINE_SNP(SNPSIFT_DBNSFP_SNP.out.vcf, 'SNP')
      SNPSIFT_EXTRACTFIELDS_SNP(CAT_ONEPERLINE_SNP.out.vcf)
    // INDEL
      COSMIC_ANNOTATION_INDEL(CAT_ANNOTATE_SNP.out.vcf)
      SNPEFF_HUMAN_INDEL(COSMIC_ANNOTATION_SNP.out.vcf, 'INDEL')
      SNPSIFT_DBNSFP_INDEL(SNPEFF_HUMAN_SNP.out.vcf, 'INDEL')
      CAT_ONEPERLINE_INDEL(SNPSIFT_DBNSFP_SNP.out.vcf, 'INDEL')
      SNPSIFT_EXTRACTFIELDS_INDEL(CAT_ONEPERLINE_SNP.out.vcf)

  // Merge SNP and INDEL and Aggregate Stats
    GATK_MERGEVCF(CAT_ONEPERLINE_SNP.out.vcf,
                  CAT_ONEPERLINE_INDEL.out.vcf)
    AGGREGATE_STATS_HUMAN(QUALITY_STATISTICS.out.quality_stats,
                          PICARD_COLLECTALIGNMENTSUMARYMETRICS.out.txt,
                          PICARD_MARKDUPLICATES.out.dedup_metrics)

  }
  
  // If Mouse
  if (params.gen_org=='mouse'){
    SNPEFF(CAT_ANNOTATE_SNP.out.vcf)
    GATK_VARIANTANNOTATOR(CAT_ANNOTATE_SNP.out.vcf,
                          SNPEFF.out.vcf)
    SNPSIFT_EXTRACTFIELDS(GATK_VARIANTANNOTATOR.out.vcf)
    AGGREGATE_STATS_MOUSE(QUALITY_STATISTICS.out.quality_stats,
                          PICARD_COLLECTALIGNMENTSUMARYMETRICS.out.txt)
  }
}
