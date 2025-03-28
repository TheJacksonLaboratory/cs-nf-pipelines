#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/wgs.nf"
include {param_log} from "${projectDir}/bin/log/wgs.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv.nf"
include {FILE_DOWNLOAD} from "${projectDir}/subworkflows/aria_download_parse"
include {CONCATENATE_LOCAL_FILES} from "${projectDir}/subworkflows/concatenate_local_files"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"

include {CLUMPIFY} from "${projectDir}/modules/bbmap/bbmap_clumpify"
include {FASTP} from "${projectDir}/modules/fastp/fastp"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"

include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {BWA_MEM_HLA} from "${projectDir}/modules/bwa/bwa_mem_hla"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {SAMTOOLS_MERGE} from "${projectDir}/modules/samtools/samtools_merge"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"

include {GATK_BASERECALIBRATOR} from "${projectDir}/modules/gatk/gatk_baserecalibrator_interval"
include {GATK_GATHERBQSRREPORTS} from "${projectDir}/modules/gatk/gatk_gatherbqsrreports"
include {GATK_APPLYBQSR} from "${projectDir}/modules/gatk/gatk_applybqsr"

include {JVARKIT_COVERAGE_CAP} from "${projectDir}/modules/jvarkit/jvarkit_biostar154220"
include {SAMTOOLS_INDEX} from "${projectDir}/modules/samtools/samtools_index"

include {PICARD_COLLECTALIGNMENTSUMMARYMETRICS} from "${projectDir}/modules/picard/picard_collectalignmentsummarymetrics"
include {PICARD_COLLECTWGSMETRICS} from "${projectDir}/modules/picard/picard_collectwgsmetrics"

include {GATK_HAPLOTYPECALLER_INTERVAL;
         GATK_HAPLOTYPECALLER_INTERVAL as GATK_HAPLOTYPECALLER_INTERVAL_GVCF} from "${projectDir}/modules/gatk/gatk_haplotypecaller_interval"
include {MAKE_VCF_LIST} from "${projectDir}/modules/utility_modules/make_vcf_list"
include {GATK_MERGEVCF_LIST} from "${projectDir}/modules/gatk/gatk_mergevcf_list"
include {GATK_COMBINEGVCFS} from "${projectDir}/modules/gatk/gatk_combinegvcfs"

include {GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from "${projectDir}/modules/gatk/gatk_selectvariants"
include {GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from "${projectDir}/modules/gatk/gatk_variantfiltration"
include {GATK_MERGEVCF;
         GATK_MERGEVCF as GATK_MERGEVCF_UNANNOTATED;
         GATK_MERGEVCF as GATK_MERGEVCF_ANNOTATED} from "${projectDir}/modules/gatk/gatk_mergevcf"

include {SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_COSMIC;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_SNP_DBSNP;
         SNPSIFT_ANNOTATE as SNPSIFT_ANNOTATE_INDEL_DBSNP} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {SNPEFF;
         SNPEFF as SNPEFF_SNP;
         SNPEFF as SNPEFF_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_snpeff"
include {SNPEFF_ONEPERLINE;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_SNP;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_oneperline"
include {SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpsift_dbnsfp"
include {SNPSIFT_EXTRACTFIELDS} from "${projectDir}/modules/snpeff_snpsift/snpsift_extractfields"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

if (params.download_data && !params.csv_input) {
    exit 1, "Data download was specified with `--download_data`. However, no input CSV file was specified with `--csv_input`. This is an invalid parameter combination. `--download_data` requires a CSV manifest. See `--help` for information."
}

// prepare reads channel
if (params.csv_input) {

    ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))

    if (params.read_type == 'PE'){
        ch_input_sample.map{it -> [it[0], [it[2], it[3]]]}.set{read_ch}
        ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
    } else if (params.read_type == 'SE') {
        ch_input_sample.map{it -> [it[0], it[2]]}.set{read_ch}
        ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

} else if (params.concat_lanes){
  
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
    // if channel is empty give error message and exit
    read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern} and file extension: ${params.extension}"}

} else {
  
  if (params.read_type == 'PE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
  }
    // if channel is empty give error message and exit
    read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern} and file extension: ${params.extension}"}
}

// main workflow
workflow WGS {
  // Step 0: Download data and concat Fastq files if needed. 
  if (params.download_data){
      FILE_DOWNLOAD(ch_input_sample)

      FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[2]]}.set{read_ch}
      FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
  }

  // Step 00: Concat local Fastq files from CSV input if required.
  if (!params.download_data && params.csv_input){
      CONCATENATE_LOCAL_FILES(ch_input_sample)
      CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[2]]}.set{read_ch}
      CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
  }
  
  // Step 00: Concat local Fastq files if required.
  if (params.concat_lanes && !params.csv_input){
      if (params.read_type == 'PE'){
          CONCATENATE_READS_PE(read_ch)
          read_ch = CONCATENATE_READS_PE.out.concat_fastq
      } else if (params.read_type == 'SE'){
          CONCATENATE_READS_SE(read_ch)
          read_ch = CONCATENATE_READS_SE.out.concat_fastq
      }
  }

  // Optional Step -- Clumpify
  if (params.deduplicate_reads) {
      CLUMPIFY(read_ch)
      trimmer_input = CLUMPIFY.out.clumpy_fastq
  } else {
      trimmer_input = read_ch
  }

  // Step 1: Read quality and adapter trimming
  FASTP(trimmer_input)
    
  FASTQC(FASTP.out.trimmed_fastq)

  // Step 2: Get Read Group Information
  READ_GROUPS(FASTP.out.trimmed_fastq, "gatk")

  if (params.split_fastq) {
    if (params.read_type == 'PE') {
      split_fastq_files = FASTP.out.trimmed_fastq
                         .map{it -> [it[0], it[1][0], it[1][1]]}
                         .splitFastq(by: params.split_fastq_bin_size, file: true, pe: true)
                         .map{it -> [it[0], [it[1], it[2]], it[1].name.split('\\.')[-2]]}
                         .combine(READ_GROUPS.out.read_groups, by: 0)
                         // from fastp the naming convention will always be *R*.fastq. 
                         // splitFastq adds an increment between *R* and .fastq. 
                         // This can be used to set an 'index' value to make file names unique. 
    } else {
      split_fastq_files = FASTP.out.trimmed_fastq
                         .map{it -> [it[0], it[1]]}
                         .splitFastq(by: params.split_fastq_bin_size, file: true)
                         .map{it -> [it[0], it[1], it[1].name.split('\\.')[-2]]}
                         .combine(READ_GROUPS.out.read_groups, by: 0)
                         // from fastp the naming convention will always be *R*.fastq. 
                         // splitFastq adds an increment between *R* and .fastq. 
                         // This can be used to set an 'index' value to make file names unique.
    }
    split_fastq_count = split_fastq_files
                    .groupTuple()
                    .map{sample, reads, index, read_group -> [sample, groupKey(sample, index.size())]}
                  
    bwa_mem_mapping = split_fastq_count
                .combine(split_fastq_files, by:0)
                .map{it -> [it[1], it[2], it[3], it[4]] }
  } else {
    bwa_mem_mapping = FASTP.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
                      .map{it -> [it[0], it[1], 'aln', it[2]]}
  }

  // Step 3: BWA-MEM Alignment
  if (params.gen_org=='mouse' | params.gen_org=='other'){
    BWA_MEM(bwa_mem_mapping)
    PICARD_SORTSAM(BWA_MEM.out.sam, 'coordinate')
  }
  if (params.gen_org=='human'){ 
  	BWA_MEM_HLA(bwa_mem_mapping)
  	PICARD_SORTSAM(BWA_MEM_HLA.out.bam, 'coordinate')
  }

  if (params.split_fastq) {
    SAMTOOLS_MERGE(PICARD_SORTSAM.out.bam.groupTuple(), 'merged_file')
    bam_file = SAMTOOLS_MERGE.out.bam
  } else {
    bam_file = PICARD_SORTSAM.out.bam
  }

  // Step 4: Variant Preprocessing - Part 1
  PICARD_MARKDUPLICATES(bam_file)

  // If Human
  ch_GATK_BASERECALIBRATOR_multiqc = Channel.empty() //optional log file for human only.
  if (params.gen_org=='human'){

    // Read a list of contigs from parameters to provide to GATK as intervals
    chroms = Channel
        .fromPath("${params.chrom_contigs}")
        .splitText()
        .map{it -> it.trim()}
    num_chroms = file(params.chrom_contigs).countLines().toInteger()

    // Calculate BQSR, scattered by chrom. gather reports and pass to applyBQSR
    GATK_BASERECALIBRATOR(PICARD_MARKDUPLICATES.out.dedup_bam.combine(chroms))
    GATK_GATHERBQSRREPORTS(GATK_BASERECALIBRATOR.out.table.groupTuple(size: num_chroms))
    ch_GATK_BASERECALIBRATOR_multiqc = GATK_GATHERBQSRREPORTS.out.table // set log file for multiqc

    // Apply BQSR
    apply_bqsr = PICARD_MARKDUPLICATES.out.dedup_bam.join(GATK_GATHERBQSRREPORTS.out.table)
    GATK_APPLYBQSR(apply_bqsr)

    if (params.coverage_cap) {
        JVARKIT_COVERAGE_CAP(GATK_APPLYBQSR.out.bam)
        SAMTOOLS_INDEX(JVARKIT_COVERAGE_CAP.out.bam)

        bam_file = JVARKIT_COVERAGE_CAP.out.bam
        index_file = SAMTOOLS_INDEX.out.bai
    } else {
        bam_file = GATK_APPLYBQSR.out.bam
        index_file = GATK_APPLYBQSR.out.bai
    }

    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(bam_file)
    PICARD_COLLECTWGSMETRICS(bam_file)

    // HaplotypeCaller does not have multithreading so it runs faster when individual chromosomes called instead of Whole Genome
    // Applies scatter intervals from above to the BQSR bam file
    chrom_channel = bam_file.join(index_file).combine(chroms)
    
    // Use the Channel in HaplotypeCaller
    GATK_HAPLOTYPECALLER_INTERVAL(chrom_channel, '')
    // Gather intervals from scattered HaplotypeCaller operations into one
    // common stream for output

    MAKE_VCF_LIST(GATK_HAPLOTYPECALLER_INTERVAL.out.vcf.groupTuple(size: num_chroms),chroms.toList())
    GATK_MERGEVCF_LIST(MAKE_VCF_LIST.out.list)

    if (params.run_gvcf) {
      // Use the Channel in HaplotypeCaller_GVCF
      GATK_HAPLOTYPECALLER_INTERVAL_GVCF(chrom_channel,'gvcf')
      GATK_COMBINEGVCFS(GATK_HAPLOTYPECALLER_INTERVAL_GVCF.out.vcf.groupTuple(size: num_chroms))
    }
  }

  // If Mouse
  if (params.gen_org=='mouse' | params.gen_org=='other'){

    if (params.coverage_cap) {
        JVARKIT_COVERAGE_CAP(PICARD_MARKDUPLICATES.out.dedup_bam)
        SAMTOOLS_INDEX(JVARKIT_COVERAGE_CAP.out.bam)

        bam_file = JVARKIT_COVERAGE_CAP.out.bam
        index_file = SAMTOOLS_INDEX.out.bai
    } else {
        bam_file = PICARD_MARKDUPLICATES.out.dedup_bam
        index_file = PICARD_MARKDUPLICATES.out.dedup_bai
    }

    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(bam_file)
    PICARD_COLLECTWGSMETRICS(bam_file)

    // Read a list of contigs from parameters to provide to GATK as intervals
    // for HaplotypeCaller variant regions
    chroms = Channel
     .fromPath("${params.chrom_contigs}")
     .splitText()
     .map{it -> it.trim()}
    
    num_chroms = file(params.chrom_contigs).countLines().toInteger()
    // number of intervals split on during calling. A 'value' variable used in groupTuple size statement. 

    // Applies scatter intervals from above to the markdup bam file
    chrom_channel = bam_file.join(index_file).combine(chroms)

    // Use the Channel in HaplotypeCaller
    GATK_HAPLOTYPECALLER_INTERVAL(chrom_channel, '')
    // Gather intervals from scattered HaplotypeCaller operations into one
    // common stream for output
  
    MAKE_VCF_LIST(GATK_HAPLOTYPECALLER_INTERVAL.out.vcf.groupTuple(size: num_chroms), chroms.toList())
    // Sort VCF within MAKE_VCF_LIST
    GATK_MERGEVCF_LIST(MAKE_VCF_LIST.out.list)

    if (params.run_gvcf) {
      // Use the Channel in HaplotypeCaller_GVCF
      GATK_HAPLOTYPECALLER_INTERVAL_GVCF(chrom_channel,'gvcf')
      GATK_COMBINEGVCFS(GATK_HAPLOTYPECALLER_INTERVAL_GVCF.out.vcf.groupTuple(size: num_chroms))
    }
  }

  // SNP
    select_var_snp = GATK_MERGEVCF_LIST.out.vcf.join(GATK_MERGEVCF_LIST.out.idx)
    GATK_SELECTVARIANTS_SNP(select_var_snp, 'SNP', 'selected_SNP')
    var_filter_snp = GATK_SELECTVARIANTS_SNP.out.vcf.join(GATK_SELECTVARIANTS_SNP.out.idx)
    GATK_VARIANTFILTRATION_SNP(var_filter_snp, 'SNP')

  // INDEL
    select_var_indel = GATK_MERGEVCF_LIST.out.vcf.join(GATK_MERGEVCF_LIST.out.idx)
    GATK_SELECTVARIANTS_INDEL(select_var_indel, 'INDEL', 'selected_INDEL')
    var_filter_indel = GATK_SELECTVARIANTS_INDEL.out.vcf.join(GATK_SELECTVARIANTS_INDEL.out.idx)
    GATK_VARIANTFILTRATION_INDEL(var_filter_indel, 'INDEL')

  // For other genome, expectation is that dbSNP will not exist.  
  if (params.gen_org=='mouse' | params.gen_org=='human'){
    SNPSIFT_ANNOTATE_SNP_DBSNP(GATK_VARIANTFILTRATION_SNP.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')
    SNPSIFT_ANNOTATE_INDEL_DBSNP(GATK_VARIANTFILTRATION_INDEL.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')
  }

  // If Human
  if (params.gen_org=='human'){

    // SNP
      SNPSIFT_ANNOTATE_SNP_COSMIC(SNPSIFT_ANNOTATE_SNP_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmicID')
      SNPEFF_SNP(SNPSIFT_ANNOTATE_SNP_COSMIC.out.vcf, 'SNP', 'vcf')
      SNPSIFT_DBNSFP_SNP(SNPEFF_SNP.out.vcf, 'SNP')
      SNPEFF_ONEPERLINE_SNP(SNPSIFT_DBNSFP_SNP.out.vcf, 'SNP')
    // INDEL
      SNPSIFT_ANNOTATE_INDEL_COSMIC(SNPSIFT_ANNOTATE_INDEL_DBSNP.out.vcf, params.cosmic, params.cosmic_index, 'cosmicID')
      SNPEFF_INDEL(SNPSIFT_ANNOTATE_INDEL_COSMIC.out.vcf, 'INDEL', 'vcf')
      SNPSIFT_DBNSFP_INDEL(SNPEFF_INDEL.out.vcf, 'INDEL')
      SNPEFF_ONEPERLINE_INDEL(SNPSIFT_DBNSFP_INDEL.out.vcf, 'INDEL')
      
    // Merge SNP and INDEL and Aggregate Stats
      vcf_files_unannotated = SNPSIFT_ANNOTATE_SNP_COSMIC.out.vcf.join(SNPSIFT_ANNOTATE_INDEL_COSMIC.out.vcf)
      GATK_MERGEVCF_UNANNOTATED(vcf_files_unannotated, 'SNP_INDEL_filtered_unannotated_final')

      vcf_files_annotated = SNPEFF_ONEPERLINE_SNP.out.vcf.join(SNPEFF_ONEPERLINE_INDEL.out.vcf)
      GATK_MERGEVCF_ANNOTATED(vcf_files_annotated, 'SNP_INDEL_filtered_annotated_final')
      
      SNPSIFT_EXTRACTFIELDS(GATK_MERGEVCF_ANNOTATED.out.vcf)
  }

  // If Mouse
  if (params.gen_org=='mouse'){
    // Merge SNP and INDEL

    vcf_files = SNPSIFT_ANNOTATE_SNP_DBSNP.out.vcf.join(SNPSIFT_ANNOTATE_INDEL_DBSNP.out.vcf)

    GATK_MERGEVCF(vcf_files, 'SNP_INDEL_filtered_unannotated_final')

    SNPEFF(GATK_MERGEVCF.out.vcf, 'BOTH', 'vcf')

    SNPEFF_ONEPERLINE(SNPEFF.out.vcf, 'BOTH')

    SNPSIFT_EXTRACTFIELDS(SNPEFF_ONEPERLINE.out.vcf)
  }

  if (params.gen_org=='other'){
  // For other genomes, there will likely not be SNP EFF annotations, but merge still needs to happen. 
    vcf_files = GATK_VARIANTFILTRATION_SNP.out.vcf.join(GATK_VARIANTFILTRATION_INDEL.out.vcf)

    GATK_MERGEVCF(vcf_files, 'SNP_INDEL_filtered_unannotated_final')
  }

  ch_multiqc_files = Channel.empty()
  ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.quality_json.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(ch_GATK_BASERECALIBRATOR_multiqc.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.txt.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTWGSMETRICS.out.txt.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.dedup_metrics.collect{it[1]}.ifEmpty([]))

  MULTIQC (
      ch_multiqc_files.collect()
  )

}
