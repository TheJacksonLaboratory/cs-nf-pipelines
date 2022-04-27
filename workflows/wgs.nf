#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/wgs.nf'
include {param_log} from '../bin/log/wgs.nf'
include {CONCATENATE_READS_PE} from '../modules/utility_modules/concatenate_reads_PE'
include {CONCATENATE_READS_SE} from '../modules/utility_modules/concatenate_reads_SE'
include {BWA_MEM} from '../modules/bwa/bwa_mem'
include {BWA_MEM_HLA} from '../modules/bwa/bwa_mem_hla'
include {COSMIC_ANNOTATION as COSMIC_ANNOTATION_SNP;
         COSMIC_ANNOTATION as COSMIC_ANNOTATION_INDEL} from '../modules/cosmic/cosmic_annotation'
include {VCF_ANNOTATE as VCF_ANNOTATE_SNP;
         VCF_ANNOTATE as VCF_ANNOTATE_INDEL} from '../modules/vcftools/vcf_annotate'
include {SNPEFF;
         SNPEFF as SNPEFF_SNP;
         SNPEFF as SNPEFF_INDEL} from '../modules/snpeff_snpsift/snpeff_snpeff'
include {SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_SNP;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_INDEL} from '../modules/snpeff_snpsift/snpeff_oneperline'
include {SNPSIFT_EXTRACTFIELDS} from '../modules/snpeff_snpsift/snpsift_extractfields'
include {SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from '../modules/snpeff_snpsift/snpsift_dbnsfp'
include {AGGREGATE_STATS} from '../modules/utility_modules/aggregate_stats_wgs'
include {READ_GROUPS} from '../modules/utility_modules/read_groups'
include {QUALITY_STATISTICS} from '../modules/utility_modules/quality_stats'
include {PICARD_SORTSAM} from '../modules/picard/picard_sortsam'
include {PICARD_MARKDUPLICATES} from '../modules/picard/picard_markduplicates'
include {PICARD_COLLECTALIGNMENTSUMMARYMETRICS} from '../modules/picard/picard_collectalignmentsummarymetrics'
include {PICARD_COLLECTWGSMETRICS} from '../modules/picard/picard_collectwgsmetrics'
include {GATK_REALIGNERTARGETCREATOR} from '../modules/gatk/gatk_realignertargetcreator'
include {GATK_BASERECALIBRATOR} from '../modules/gatk/gatk_baserecalibrator'
include {GATK_APPLYBQSR} from '../modules/gatk/gatk_applybqsr'
include {GATK_INDELREALIGNER} from '../modules/gatk/gatk_indelrealigner'
include {GATK_MERGEVCF} from '../modules/gatk/gatk_mergevcf'
include {GATK_MERGEVCF_LIST} from '../modules/gatk/gatk_mergevcf_list'
include {GATK_VARIANTANNOTATOR} from '../modules/gatk/gatk_variantannotator'
include {GATK_HAPLOTYPECALLER_INTERVAL} from '../modules/gatk/gatk_haplotypecaller_interval'
include {GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from '../modules/gatk/gatk_selectvariants'
include {GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from '../modules/gatk/gatk_variantfiltration'
include {MAKE_VCF_LIST} from '../modules/utility_modules/make_vcf_list'

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
workflow WGS {
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
    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(GATK_APPLYBQSR.out.bam)
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
    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(GATK_INDELREALIGNER.out.bam)
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

  // Cat Output to vcf-annotate* and add dbSNP annotations. 
    VCF_ANNOTATE_SNP(GATK_VARIANTFILTRATION_SNP.out.vcf, 'SNP')
    VCF_ANNOTATE_INDEL(GATK_VARIANTFILTRATION_INDEL.out.vcf, 'INDEL')
  
// Final Post-Processing Steps Differ for Human and Mouse

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
    // Merge SNP and INDEL
    GATK_MERGEVCF(VCF_ANNOTATE_SNP.out.vcf,
                  VCF_ANNOTATE_INDEL.out.vcf)
    SNPEFF(GATK_MERGEVCF.out.vcf, 'BOTH', 'gatk')
    GATK_VARIANTANNOTATOR(GATK_MERGEVCF.out.vcf,
                          SNPEFF.out.vcf)
    SNPSIFT_EXTRACTFIELDS(GATK_VARIANTANNOTATOR.out.vcf)

  }

  // may replace with multiqc
  AGGREGATE_STATS(QUALITY_STATISTICS.out.quality_stats,
                  PICARD_MARKDUPLICATES.out.dedup_metrics,
                  PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.txt,
                  PICARD_COLLECTWGSMETRICS.out.txt)
}
