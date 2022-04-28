#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/wes.nf'
include {param_log} from '../bin/log/wes.nf'
include {getLibraryId} from '../bin/shared/getLibraryId.nf'
include {CONCATENATE_READS_PE} from '../modules/utility_modules/concatenate_reads_PE'
include {CONCATENATE_READS_SE} from '../modules/utility_modules/concatenate_reads_SE'
include {BWA_MEM} from '../modules/bwa/bwa_mem'
include {SAMTOOLS_INDEX} from '../modules/samtools/samtools_index'
include {READ_GROUPS} from '../modules/utility_modules/read_groups'
include {QUALITY_STATISTICS} from '../modules/utility_modules/quality_stats'
include {AGGREGATE_STATS} from '../modules/utility_modules/aggregate_stats_wes'
include {COSMIC_ANNOTATION;
        COSMIC_ANNOTATION as COSMIC_ANNOTATION_SNP;
        COSMIC_ANNOTATION as COSMIC_ANNOTATION_INDEL} from '../modules/cosmic/cosmic_annotation'
include {PICARD_SORTSAM} from '../modules/picard/picard_sortsam'
include {PICARD_MARKDUPLICATES} from '../modules/picard/picard_markduplicates'
include {PICARD_COLLECTHSMETRICS} from '../modules/picard/picard_collecthsmetrics'
include {SNPEFF;
         SNPEFF as SNPEFF_SNP;
         SNPEFF as SNPEFF_INDEL} from '../modules/snpeff_snpsift/snpeff_snpeff'
include {SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_SNP;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_INDEL} from '../modules/snpeff_snpsift/snpeff_oneperline'
include {SNPSIFT_EXTRACTFIELDS} from '../modules/snpeff_snpsift/snpsift_extractfields'
include {SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from '../modules/snpeff_snpsift/snpsift_dbnsfp'
include {GATK_HAPLOTYPECALLER;
         GATK_HAPLOTYPECALLER as GATK_HAPLOTYPECALLER_GVCF} from '../modules/gatk/gatk_haplotypecaller'
include {GATK_INDEXFEATUREFILE} from '../modules/gatk/gatk_indexfeaturefile'
include {GATK_VARIANTFILTRATION;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from '../modules/gatk/gatk_variantfiltration'
include {GATK_VARIANTANNOTATOR} from '../modules/gatk/gatk_variantannotator'
include {GATK_MERGEVCF} from '../modules/gatk/gatk_mergevcf'
include {GATK_SELECTVARIANTS;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from '../modules/gatk/gatk_selectvariants'
include {GATK_BASERECALIBRATOR} from '../modules/gatk/gatk_baserecalibrator'
include {GATK_APPLYBQSR} from '../modules/gatk/gatk_applybqsr'

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// prepare reads channel
if (params.concat_lanes){
  if (params.read_type == 'PE'){
    read_ch = Channel
            .fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true, flat:true )
            .map { file, file1, file2 -> tuple(getLibraryId(file), file1, file2) }
            .groupTuple()
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}", checkExists:true, size:1 )
                .map { file, file1 -> tuple(getLibraryId(file), file1) }
                .groupTuple()
                .map{t-> [t[0], t[1].flatten()]}
  }
} else {
  if (params.read_type == 'PE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
  }
}

// if channel is empty give error message and exit
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

// main workflow
workflow WES {
  // Step 0: Concatenate Fastq files if required. 
  if (params.concat_lanes){
    if (params.read_type == 'PE'){
        CONCATENATE_READS_PE(read_ch)
        read_ch = CONCATENATE_READS_PE.out.concat_fastq
    } else if (params.read_type == 'SE'){
        CONCATENATE_READS_SE(read_ch)
        read_ch = CONCATENATE_READS_SE.out.concat_fastq
    }
  }

  // Step 1: Qual_Stat
  QUALITY_STATISTICS(read_ch)

  // Step 2: Get Read Group Information
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq, "gatk")

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
        SNPEFF_SNP(COSMIC_ANNOTATION_SNP.out.vcf, 'SNP', 'vcf')
        SNPSIFT_DBNSFP_SNP(SNPEFF_SNP.out.vcf, 'SNP')
        SNPEFF_ONEPERLINE_SNP(SNPSIFT_DBNSFP_SNP.out.vcf, 'SNP')

      // INDEL
        COSMIC_ANNOTATION_INDEL(GATK_VARIANTFILTRATION_INDEL.out.vcf)
        SNPEFF_INDEL(COSMIC_ANNOTATION_INDEL.out.vcf, 'INDEL', 'vcf')
        SNPSIFT_DBNSFP_INDEL(SNPEFF_INDEL.out.vcf, 'INDEL')
        SNPEFF_ONEPERLINE_INDEL(SNPSIFT_DBNSFP_INDEL.out.vcf, 'INDEL')

    // Step 10: Post Variant Calling Processing - Part 2
      GATK_MERGEVCF(SNPEFF_ONEPERLINE_SNP.out.vcf,
                    SNPEFF_ONEPERLINE_INDEL.out.vcf)
      
      SNPSIFT_EXTRACTFIELDS(GATK_MERGEVCF.out.vcf)


  } else if (params.gen_org=='mouse'){

    // Step 6: Variant Pre-Processing - Part 3
      PICARD_COLLECTHSMETRICS(PICARD_MARKDUPLICATES.out.dedup_bam,
                              PICARD_MARKDUPLICATES.out.dedup_bai)

    // Step 7: Variant Calling
      GATK_HAPLOTYPECALLER(PICARD_MARKDUPLICATES.out.dedup_bam,
                           PICARD_MARKDUPLICATES.out.dedup_bai,
                          'variant')

    // Step 8: Variant Filtration
      GATK_VARIANTFILTRATION(GATK_HAPLOTYPECALLER.out.vcf,
                             GATK_HAPLOTYPECALLER.out.idx,
                            'BOTH')

    // Step 9: Post Variant Calling Processing
      SNPEFF(GATK_VARIANTFILTRATION.out.vcf, 'BOTH', 'gatk')
      GATK_VARIANTANNOTATOR(GATK_VARIANTFILTRATION.out.vcf,
                            SNPEFF.out.vcf)
      SNPSIFT_EXTRACTFIELDS(GATK_VARIANTANNOTATOR.out.vcf)

  }

  // Step 11: Aggregate Stats
  AGGREGATE_STATS(QUALITY_STATISTICS.out.quality_stats,
						      PICARD_COLLECTHSMETRICS.out.hsmetrics,
						      PICARD_MARKDUPLICATES.out.dedup_metrics)
}
