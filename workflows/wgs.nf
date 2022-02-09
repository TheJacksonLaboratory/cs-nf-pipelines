#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/wes.nf'
include {param_log} from '../bin/log/wes.nf'
include {BWA_MEM} from '../modules/bwa'
include {READ_GROUPS} from '../modules/read_groups'
include {QUALITY_STATISTICS} from '../modules/quality_stats'

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

  /* If Human: Step 5-X
  if (params.gen_org=='human'){
    GATK_REALIGNERTARGETCREATOR(PICARD_MARKDUPLICATES.out.bam) //may need bai
    GATK_INDELREALIGNER(GATK_REALIGNERTARGETCREATOR.out.intervals)
    BaseRecalibrator
    PrintReads
    CollectAlignmentSummaryMetrics (picard)


    HaplotypeCaller NOTE ib for the chromosomes (1-23 human dif for mouse?)
  }
  */

  else if (params.gen_org=='mouse'){
  // Step 5:
    GATK_REALIGNERTARGETCREATOR(PICARD_MARKDUPLICATES.out.bam) //may need bai
    GATK_INDELREALIGNER(GATK_REALIGNERTARGETCREATOR.out.intervals)
    PICARD_COLLECTALIGNMENTSUMARYMETRICS(GATK_INDELREALIGNER.out.bam)

    // Next Week
    // HaplotypeCaller NOTE ib for the chromosomes (1-23 human dif for mouse?)
  }

  // Mouse and Human Same Here

  //Next Week with HaplotypeCaller (will want a index output added as well)
  GATK_MERGEVCF Note: change this to take in a list (will need to update wes)

  // SNP
    GATK_SELECTVARIANTS_SNP(GATK_MERGEVCF.out.vcf,
                          GATK_MERGEVCF.out.idx,
                         'SNP')
    BCF_SORT_SNP(GATK_SELECTVARIANTS_SNP.out.vcf) // DOES THIS GIVE A INDEX? IS THIS NECESSARY?
    GATK_INDEXFEATUREFILE_SNP(BCF_SORT_SNP.out.vcf) // NECESSARY?
    GATK_VARIANTFILTRATION_SNP(BCF_SORT_SNP.out.vcf,
                               GATK_INDEXFEATUREFILE_SNP.out.idx,
                              'SNP')
  // INDEL
  GATK_SELECTVARIANTS_INDEL(GATK_MERGEVCF.out.vcf,
                        GATK_MERGEVCF.out.idx,
                       'INDEL')
  BCF_SORT_INDEL(GATK_SELECTVARIANTS_INDEL.out.vcf) // DOES THIS GIVE A INDEX? IS THIS NECESSARY?
  GATK_INDEXFEATUREFILE_INDEL(BCF_SORT_INDEL.out.vcf) // NECESSARY?
  GATK_VARIANTFILTRATION_INDEL(BCF_SORT_INDEL.out.vcf,
                             GATK_INDEXFEATUREFILE_INDEL.out.idx,
                            'INDEL')

  /* Finishing Steps Dif Mouse and Human
  if (params.gen_org=='human'){
    "Step 9: Post Variant Calling Processing, Part 2 lots to break down here"
    CombineVariants
    AGGREGATE_STATS_HUMAN
  }
  */
  if (params.gen_org=='mouse'){
    // NEXT WEEK
    "Step 9: Post Variant Calling Processing, Part 2 lots to break down here"

    SNPEFF(.out.vcf)
    GATK_VARIANTANNOTATOR(SAMPLE.out.vcf,
                          SNPEFF.out.vcf)
    SNPSIFT_EXTRACTFIELDS(GATK_VARIANTANNOTATOR.out.vcf)
    AGGREGATE_STATS_MOUSE(QUALITY_STATISTICS.out.quality_stats,
                          PICARD_COLLECTALIGNMENTSUMARYMETRICS.out.txt)
  }
}
