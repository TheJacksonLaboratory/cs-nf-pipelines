#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/wgs.nf'
include {param_log} from '../bin/log/wgs.nf'
include {BWA_MEM;
         BWA_MEM_HLA} from '../modules/bwa'
include {COSMIC_ANNOTATION as COSMIC_ANNOTATION_SNP;
         COSMIC_ANNOTATION as COSMIC_ANNOTATION_INDEL} from '../modules/cosmic'
include {VCF_ANNOTATE as VCF_ANNOTATE_SNP;
         VCF_ANNOTATE as VCF_ANNOTATE_INDEL} from '../bin/wgs/vcf_annotate'
include {SNPEFF;
         SNPEFF as SNPEFF_SNP;
         SNPEFF as SNPEFF_INDEL;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_SNP;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_INDEL} from '../modules/snpeff'
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
         PICARD_COLLECTALIGNMENTSUMARYMETRICS;
         PICARD_COLLECTWGSMETRICS} from '../modules/picard'
include {GATK_REALIGNERTARGETCREATOR;
         GATK_BASERECALIBRATOR;
         GATK_APPLYBQSR;
         GATK_INDELREALIGNER;
         GATK_MERGEVCF;
         GATK_MERGEVCF_LIST;
         GATK_VARIANTANNOTATOR;
         GATK_HAPLOTYPECALLER_INTERVAL;
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
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq, "gatk")
  // Step 3: BWA-MEM Alignment
  if (params.gen_org=='mouse'){
    BWA_MEM(QUALITY_STATISTICS.out.trimmed_fastq, READ_GROUPS.out.read_groups)
    PICARD_SORTSAM(BWA_MEM.out.sam)
  }
  if (params.gen_org=='human'){ 
  	BWA_MEM_HLA(QUALITY_STATISTICS.out.trimmed_fastq, READ_GROUPS.out.read_groups)
  	PICARD_SORTSAM(BWA_MEM_HLA.out.bam)
  }
  // Step 4: Variant Preprocessing - Part 1
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
  // Step 5 Depricated in GATK 4
  GATK_REALIGNERTARGETCREATOR(PICARD_MARKDUPLICATES.out.dedup_bam)
  GATK_INDELREALIGNER(PICARD_MARKDUPLICATES.out.dedup_bam,
                      GATK_REALIGNERTARGETCREATOR.out.intervals)
  // If Human
  if (params.gen_org=='human'){
    GATK_BASERECALIBRATOR(GATK_INDELREALIGNER.out.bam)
    GATK_APPLYBQSR(GATK_INDELREALIGNER.out.bam,
                    GATK_BASERECALIBRATOR.out.table)
    PICARD_COLLECTALIGNMENTSUMARYMETRICS(GATK_APPLYBQSR.out.bam)
    PICARD_COLLECTWGSMETRICS(GATK_APPLYBQSR.out.bam)

    // Create a chromosome channel. HaplotypeCaller does not have multithreading so it runs faster when individual chromosomes called instead of Whole Genome
    data = GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai)
    
    // Read a list of contigs from parameters to provide to GATK as intervals
    // for HaplotypeCaller variant regions
    chroms = Channel
     .fromPath("${params.chrom_contigs}")
     .splitText()
     .map{it -> it.trim()}
    
    // Applies scatter intervals from above to the BQSR bam file
    chrom_channel = data.combine(chroms)
    
    // Use the Channel in HaplotypeCaller
    GATK_HAPLOTYPECALLER_INTERVAL(chrom_channel)
    // Gather intervals from scattered HaplotypeCaller operations into one
    // common stream for output
    MAKE_VCF_LIST(GATK_HAPLOTYPECALLER_INTERVAL.out.vcf.groupTuple(),chroms.toList())
    GATK_MERGEVCF_LIST(MAKE_VCF_LIST.out.list)
  }

  // If Mouse
  if (params.gen_org=='mouse'){
    PICARD_COLLECTALIGNMENTSUMARYMETRICS(GATK_INDELREALIGNER.out.bam)
    PICARD_COLLECTWGSMETRICS(GATK_INDELREALIGNER.out.bam)

    // create a chromosome channel. HaplotypeCaller runs faster when individual chromosomes called instead of Whole Genome
    data = GATK_INDELREALIGNER.out.bam.join(GATK_INDELREALIGNER.out.bai)
    
    // Read a list of contigs from parameters to provide to GATK as intervals
    // for HaplotypeCaller variant regions
    chroms = Channel
     .fromPath("${params.chrom_contigs}")
     .splitText()
     .map{it -> it.trim()}
    
    // Applies scatter intervals from above to the BQSR bam file
    chrom_channel = data.combine(chroms)

    // Use the Channel in HaplotypeCaller
    GATK_HAPLOTYPECALLER_INTERVAL(chrom_channel)
    // Gather intervals from scattered HaplotypeCaller operations into one
    // common stream for output
    MAKE_VCF_LIST(GATK_HAPLOTYPECALLER_INTERVAL.out.vcf.groupTuple(), chroms.toList())
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
      SNPEFF_SNP(COSMIC_ANNOTATION_SNP.out.vcf, 'SNP', 'vcf')
      SNPSIFT_DBNSFP_SNP(SNPEFF_SNP.out.vcf, 'SNP')
      SNPEFF_ONEPERLINE_SNP(SNPSIFT_DBNSFP_SNP.out.vcf, 'SNP')
    // INDEL
      COSMIC_ANNOTATION_INDEL(VCF_ANNOTATE_INDEL.out.vcf)
      SNPEFF_INDEL(COSMIC_ANNOTATION_INDEL.out.vcf, 'INDEL', 'vcf')
      SNPSIFT_DBNSFP_INDEL(SNPEFF_INDEL.out.vcf, 'INDEL')
      SNPEFF_ONEPERLINE_INDEL(SNPSIFT_DBNSFP_INDEL.out.vcf, 'INDEL')
      
  // Merge SNP and INDEL and Aggregate Stats
    GATK_MERGEVCF(SNPEFF_ONEPERLINE_SNP.out.vcf,
                  SNPEFF_ONEPERLINE_INDEL.out.vcf)

    SNPSIFT_EXTRACTFIELDS(GATK_MERGEVCF.out.vcf)
  }

  // If Mouse
  if (params.gen_org=='mouse'){
    SNPEFF(VCF_ANNOTATE_SNP.out.vcf, 'BOTH', 'gatk')
    GATK_VARIANTANNOTATOR(VCF_ANNOTATE_SNP.out.vcf,
                          SNPEFF.out.vcf)
    SNPSIFT_EXTRACTFIELDS(GATK_VARIANTANNOTATOR.out.vcf)

  }

  // may replace with multiqc
  AGGREGATE_STATS(QUALITY_STATISTICS.out.quality_stats,
                  PICARD_MARKDUPLICATES.out.dedup_metrics,
                  PICARD_COLLECTALIGNMENTSUMARYMETRICS.out.txt,
                  PICARD_COLLECTWGSMETRICS.out.txt)
}
