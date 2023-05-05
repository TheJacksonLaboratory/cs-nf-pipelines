#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/amplicon.nf"
include {param_log} from "${projectDir}/bin/log/amplicon.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {TRIM_FASTQ as CUTADAPT} from "${projectDir}/modules/cutadapt/cutadapt_trim_fastq"
include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {SAMTOOLS_SORT as SAMTOOLS_SORT_PRIMERCLIP;
         SAMTOOLS_SORT as SAMTOOLS_SORT_CALLING} from "${projectDir}/modules/samtools/samtools_sort"
include {PRIMERCLIP} from "${projectDir}/modules/primerclip/primerclip"
include {TARGET_COVERAGE_METRICS} from "${projectDir}/modules/bedtools/bedtools_amplicon_metrics"
include {SNPSIFT_ANNOTATE} from "${projectDir}/modules/snpeff_snpsift/snpsift_annotate"
include {PICARD_COLLECTTARGETPCRMETRICS} from "${projectDir}/modules/picard/picard_collecttargetpcrmetrics"
include {GATK_BASERECALIBRATOR} from "${projectDir}/modules/gatk/gatk_baserecalibrator"
include {GATK_APPLYBQSR} from "${projectDir}/modules/gatk/gatk_applybqsr"
include {GATK_HAPLOTYPECALLER} from "${projectDir}/modules/gatk/gatk_haplotypecaller"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"

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

workflow AMPLICON {
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

    CUTADAPT(read_ch)

    FASTQC(CUTADAPT.out.paired_trimmed_fastq)

    // Step 2: Get Read Group Information
    READ_GROUPS(CUTADAPT.out.paired_trimmed_fastq, "gatk")

    // Step 3: BWA-MEM Alignment
    bwa_mem_mapping = CUTADAPT.out.paired_trimmed_fastq.join(READ_GROUPS.out.read_groups)
    BWA_MEM(bwa_mem_mapping)

    SAMTOOLS_SORT_PRIMERCLIP(BWA_MEM.out.sam, '-O sam -n', 'sam')

    PRIMERCLIP(SAMTOOLS_SORT_PRIMERCLIP.out.sorted_file)

    SAMTOOLS_SORT_CALLING(PRIMERCLIP.out.sam, '-O bam', 'bam')

    PICARD_COLLECTTARGETPCRMETRICS(SAMTOOLS_SORT_CALLING.out.sorted_file)

    TARGET_COVERAGE_METRICS(SAMTOOLS_SORT_CALLING.out.sorted_file)

    /*
    Important: While the use of the Picard tool, MarkDuplicates, is a common quality control step to identify
    low-complexity libraries, MarkDuplicates cannot be used on data derived from PCR-based target enrichment
    methods such as the xGen Amplicon Panels. Since these targeted panels contain high numbers of identical
    library fragments (particularly regarding alignment start position), MarkDuplicates cannot appropriately 
    analyze Amplicon libraries.
    https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/application-note/primerclip-a-tool-for-trimming-primer-sequences-application-note.pdf?sfvrsn=cf83e107_14
    */

    GATK_BASERECALIBRATOR(SAMTOOLS_SORT_CALLING.out.sorted_file)
    
    GATK_APPLYBQSR(SAMTOOLS_SORT_CALLING.out.sorted_file.join(GATK_BASERECALIBRATOR.out.table))

    GATK_HAPLOTYPECALLER(GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai), '')

    SNPSIFT_ANNOTATE(GATK_HAPLOTYPECALLER.out.vcf, params.dbSNP, params.dbSNP_index, 'dbsnpID')

    // MultiQC
    // coverage metrics? 
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.cutadapt_log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GATK_BASERECALIBRATOR.out.table.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTTARGETPCRMETRICS.out.txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PRIMERCLIP.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TARGET_COVERAGE_METRICS.out.qc_metrics.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
}