#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/wgs.nf'
include {param_log} from '../bin/log/wgs.nf'
include {BWA_MEM} from '../modules/bwa'
include {COSMIC_ANNOTATION as COSMIC_ANNOTATION_SNP;
         COSMIC_ANNOTATION as COSMIC_ANNOTATION_INDEL} from '../modules/cosmic'
include {SNPEFF_HUMAN as SNPEFF_HUMAN_SNP;
         SNPEFF_HUMAN as SNPEFF_HUMAN_INDEL} from '../modules/snpeff'
include {VCF_ANNOTATE as VCF_ANNOTATE_SNP;
         VCF_ANNOTATE as VCF_ANNOTATE_INDEL;
         CAT_ONEPERLINE as CAT_ONEPERLINE_SNP;
         CAT_ONEPERLINE as CAT_ONEPERLINE_INDEL} from '../bin/wgs/vcf_annotate'
include {SNPEFF} from '../modules/snpeff'
include {SNPSIFT_EXTRACTFIELDS;
         SNPSIFT_EXTRACTFIELDS as SNPSIFT_EXTRACTFIELDS_SNP;
         SNPSIFT_EXTRACTFIELDS as SNPSIFT_EXTRACTFIELDS_INDEL;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from '../modules/snpsift'
include {AGGREGATE_STATS} from '../bin/wgs/aggregate_stats_wgs'
include {READ_GROUPS} from '../modules/read_groups'
include {QUALITY_STATISTICS} from '../modules/quality_stats'
include {PICARD_SORTSAM;
         PICARD_MARKDUPLICATES;
         PICARD_COLLECTALIGNMENTSUMARYMETRICS} from '../modules/picard'
include {GATK_REALIGNERTARGETCREATOR;
         GATK_BASERECALIBRATOR;
         GATK_APPLYBQSR;
         GATK_INDELREALIGNER;
         GATK_MERGEVCF;
         GATK_MERGEVCF_LIST;
         GATK_VARIANTANNOTATOR;
         GATK_HAPLOTYPECALLER_WGS;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from '../modules/gatk'
include {MAKE_VCF_LIST} from '../bin/wgs/make_vcf_list'

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
                      GATK_REALIGNERTARGETCREATOR.out.intervals) // always saved mouse optional for human (like markdups)
  // If Human
  if (params.gen_org=='human'){
    GATK_BASERECALIBRATOR(GATK_INDELREALIGNER.out.bam)
    GATK_APPLYBQSR(GATK_INDELREALIGNER.out.bam,
                    GATK_BASERECALIBRATOR.out.table)
    PICARD_COLLECTALIGNMENTSUMARYMETRICS(GATK_APPLYBQSR.out.bam)

    // create a chromosome channel. HaplotypeCaller runs faster when individual chromosomes called instead of Whole Genome
    data = GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai)
    chromes = Channel.of('chr1', 'chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                         'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
                         'chr21','chr22','chrM','chrX','chrY')
    chrome_channel = data.combine(chromes)

    // Use the Channel in HaplotypeCaller
    GATK_HAPLOTYPECALLER_WGS(chrome_channel) // GATK_HAPLOTYPECALLER_INTERVAL
    MAKE_VCF_LIST(GATK_HAPLOTYPECALLER_WGS.out.vcf.groupTuple())
    GATK_MERGEVCF_LIST(MAKE_VCF_LIST.out.list)
  }

  // If Mouse
  if (params.gen_org=='mouse'){
    PICARD_COLLECTALIGNMENTSUMARYMETRICS(GATK_INDELREALIGNER.out.bam)

    // create a chromosome channel. HaplotypeCaller runs faster when individual chromosomes called instead of Whole Genome
    data = GATK_INDELREALIGNER.out.bam.join(GATK_INDELREALIGNER.out.bai)
    chromes = Channel.of('chr1', 'chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                         'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
                         'chrM','chrX','chrY')
    chrome_channel = data.combine(chromes)

    // Use the Channel in HaplotypeCaller
    GATK_HAPLOTYPECALLER_WGS(chrome_channel)
    MAKE_VCF_LIST(GATK_HAPLOTYPECALLER_WGS.out.vcf.groupTuple())
    // Sort VCF within MAKE_VCF_LIST
    GATK_MERGEVCF_LIST(MAKE_VCF_LIST.out.list)
    }

  // SNP
    GATK_SELECTVARIANTS_SNP(GATK_MERGEVCF_LIST.out.vcf,
                            GATK_MERGEVCF_LIST.out.idx,
                           'SNP')
    GATK_VARIANTFILTRATION_SNP(GATK_SELECTVARIANTS_SNP.out.vcf,
                               GATK_SELECTVARIANTS_SNP.out.idx,
                              'SNP')
  // INDEL
    GATK_SELECTVARIANTS_INDEL(GATK_MERGEVCF_LIST.out.vcf,
                              GATK_MERGEVCF_LIST.out.idx,
                             'INDEL')
    GATK_VARIANTFILTRATION_INDEL(GATK_SELECTVARIANTS_INDEL.out.vcf,
                                 GATK_SELECTVARIANTS_INDEL.out.idx,
                                'INDEL')

  // Cat Output to vcf-annotate*
    VCF_ANNOTATE_SNP(GATK_VARIANTFILTRATION_SNP.out.vcf)
    VCF_ANNOTATE_INDEL(GATK_VARIANTFILTRATION_INDEL.out.vcf)

// Final Post-Processing Steps Slightly Different for Mouse and Human

  // If Human
  if (params.gen_org=='human'){

    // SNP
      COSMIC_ANNOTATION_SNP(VCF_ANNOTATE_SNP.out.vcf)
      SNPEFF_HUMAN_SNP(COSMIC_ANNOTATION_SNP.out.vcf, 'SNP')
      SNPSIFT_DBNSFP_SNP(SNPEFF_HUMAN_SNP.out.vcf, 'SNP')
      CAT_ONEPERLINE_SNP(SNPSIFT_DBNSFP_SNP.out.vcf, 'SNP')
      SNPSIFT_EXTRACTFIELDS_SNP(CAT_ONEPERLINE_SNP.out.vcf)
    // INDEL
      COSMIC_ANNOTATION_INDEL(VCF_ANNOTATE_INDEL.out.vcf)
      SNPEFF_HUMAN_INDEL(COSMIC_ANNOTATION_INDEL.out.vcf, 'INDEL')
      SNPSIFT_DBNSFP_INDEL(SNPEFF_HUMAN_INDEL.out.vcf, 'INDEL')
      CAT_ONEPERLINE_INDEL(SNPSIFT_DBNSFP_INDEL.out.vcf, 'INDEL')
      SNPSIFT_EXTRACTFIELDS_INDEL(CAT_ONEPERLINE_INDEL.out.vcf)

  // Merge SNP and INDEL and Aggregate Stats
    GATK_MERGEVCF(CAT_ONEPERLINE_SNP.out.vcf,
                  CAT_ONEPERLINE_INDEL.out.vcf)
  }

  // If Mouse
  if (params.gen_org=='mouse'){
    SNPEFF(VCF_ANNOTATE_SNP.out.vcf)
    GATK_VARIANTANNOTATOR(VCF_ANNOTATE_SNP.out.vcf,
                          SNPEFF.out.vcf)
    SNPSIFT_EXTRACTFIELDS(GATK_VARIANTANNOTATOR.out.vcf)

  }

  // may replace with multiqc
  AGGREGATE_STATS(QUALITY_STATISTICS.out.quality_stats,
                  PICARD_MARKDUPLICATES.out.dedup_metrics,
                  PICARD_COLLECTALIGNMENTSUMARYMETRICS.out.txt)
}
