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
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {PRIMERCLIP} from "${projectDir}/modules/primerclip/primerclip"
include {PICARD_ADDORREPLACEREADGROUPS} from "${projectDir}/modules/picard/picard_addorreplacereadgroups"
include {GATK_BASERECALIBRATOR} from "${projectDir}/modules/gatk/gatk_baserecalibrator"
include {GATK_APPLYBQSR} from "${projectDir}/modules/gatk/gatk_applybqsr"
include {GATK_HAPLOTYPECALLER} from "${projectDir}/modules/gatk/gatk_haplotypecaller"


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

    PRIMERCLIP(BWA_MEM.out.sam)

    PRIMERCLIP.out.sam.join(READ_GROUPS.out.read_groups)

    PICARD_ADDORREPLACEREADGROUPS(PRIMERCLIP.out.sam) // is this step needed as read group is used in BWA? 

    // Picard CollectTargetedPcrMetrics : this module needs to be made and figured out https://gatk.broadinstitute.org/hc/en-us/articles/360037438131-CollectTargetedPcrMetrics-Picard-

    // Bedtools coverageBed : this module needs to be made, and figured out what to do with the output. 
    /*
    Important: While the use of the Picard tool, MarkDuplicates, is a common quality control step to identify
    low-complexity libraries, MarkDuplicates cannot be used on data derived from PCR-based target enrichment
    methods such as the xGen Amplicon Panels. Since these targeted panels contain high numbers of identical
    library fragments (particularly regarding alignment start position), MarkDuplicates cannot appropriately 
    analyze Amplicon libraries.
    https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/application-note/primerclip-a-tool-for-trimming-primer-sequences-application-note.pdf?sfvrsn=cf83e107_14
    */

    GATK_BASERECALIBRATOR(PICARD_ADDORREPLACEREADGROUPS.out.bam.join(PICARD_ADDORREPLACEREADGROUPS.out.bai))
    
    GATK_APPLYBQSR(PICARD_ADDORREPLACEREADGROUPS.out.bam.join(GATK_BASERECALIBRATOR.out.table))

    GATK_HAPLOTYPECALLER(GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai))

    // multiqc input
    // FASTQC.out.quality_stats
    // CUTADAPT.out.cutadapt_log
    // CollectTargetedPcrMetrics
    // coverage metrics? 
    // primer clip likely has an output log. Once we can run, this can be collected. However, primerclip is not in multiQC.

}