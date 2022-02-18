#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/wes.nf'
include {param_log} from '../bin/log/wes.nf'
include {BWA_MEM} from '../modules/bwa'
include {CAT_ANNOTATE as CAT_ANNOTATE_SNP;
         CAT_ANNOTATE as CAT_ANNOTATE_INDEL} from '../bin/wgs/cat'
include {SNPEFF} from '../modules/snpeff'
include {SNPSIFT_EXTRACTFIELDS} from '../modules/snpsift'
include {AGGREGATE_STATS_MOUSE} from '../bin/wes/aggregate_stats'
include {READ_GROUPS} from '../modules/read_groups'
include {QUALITY_STATISTICS} from '../modules/quality_stats'
include {PICARD_SORTSAM;
         PICARD_MARKDUPLICATES;
         PICARD_COLLECTALIGNMENTSUMARYMETRICS} from '../modules/picard'
include {GATK_REALIGNERTARGETCREATOR;
         GATK_BASERECALIBRATOR;
         GATK_PRINTREADS;
         GATK_INDELREALIGNER;
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
    
//    GATK_BASERECALIBRATOR(GATK_INDELREALIGNER.out.bam)
// ran into major issue here
//    GATK_PRINTREADS(GATK_INDELREALIGNER.out.bam,
//                    GATK_BASERECALIBRATOR.out.table)
  //for now skipping printreads and using indelrealigner output instead
    PICARD_COLLECTALIGNMENTSUMARYMETRICS(GATK_INDELREALIGNER.out.bam)
    GATK_HAPLOTYPECALLER_WGS(GATK_INDELREALIGNER.out.bam,
                             GATK_INDELREALIGNER.out.bai)
  }


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

/*   Finishing Steps Dif Mouse and Human
  if (params.gen_org=='human'){
    "Step 9: Post Variant Calling Processing, Part 2 lots to break down here"
    CombineVariants
    AGGREGATE_STATS_HUMAN
  }
  
  if (params.gen_org=='mouse'){
    CAT_ANNOTATE_SNP(GATK_VARIANTFILTRATION_SNP.out.vcf)
    CAT_ANNOTATE_INDEL(GATK_VARIANTFILTRATION_INDEL.out.vcf)
    SNPEFF(CAT_ANNOTATE_SNP.out.vcf)
    GATK_VARIANTANNOTATOR(CAT_ANNOTATE_SNP.out.vcf,
                          SNPEFF.out.vcf)
    SNPSIFT_EXTRACTFIELDS(GATK_VARIANTANNOTATOR.out.vcf)
    AGGREGATE_STATS_MOUSE(QUALITY_STATISTICS.out.quality_stats,
                          PICARD_COLLECTALIGNMENTSUMARYMETRICS.out.txt)
  }
*/
}
