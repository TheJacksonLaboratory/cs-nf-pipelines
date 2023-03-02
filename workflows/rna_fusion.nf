#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/rna_fusion.nf"
include {param_log} from "${projectDir}/bin/log/rna_fusion.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv.nf"
include {FILE_DOWNLOAD} from "${projectDir}/subworkflows/aria_download_parse"
include {CONCATENATE_LOCAL_FILES} from "${projectDir}/subworkflows/concatenate_local_files"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {GUNZIP} from "${projectDir}/modules/utility_modules/gunzip"
include {XENOME_CLASSIFY} from "${projectDir}/modules/xenome/xenome"
include {FASTQ_PAIR} from "${projectDir}/modules/fastq-tools/fastq-pair"
include {FASTQ_SORT as FASTQ_SORT_HUMAN;
         FASTQ_SORT as FASTQ_SORT_MOUSE} from "${projectDir}/modules/fastq-tools/fastq-sort"

include {STAR_ALIGN as STAR_ARRIBA;
         STAR_ALIGN as STAR_SQUID} from "${projectDir}/modules/star/star_align"

include {SAMTOOLS_SORT as SORT_ARRIBA;
         SAMTOOLS_SORT as SORT_SQUID} from "${projectDir}/modules/samtools/samtools_sort_only"
include {SAMTOOLS_INDEX as INDEX_ARRIBA} from "${projectDir}/modules/samtools/samtools_index"

include {ARRIBA} from "${projectDir}/modules/arriba/arriba"

include {FUSIONCATCHER} from "${projectDir}/modules/fusioncatcher/fusioncatcher"

include {JAFFA} from "${projectDir}/modules/jaffa/jaffa"

include {KALLISTO_QUANT} from "${projectDir}/modules/kallisto/kallisto_quant"
include {KALLISTO_INSERT_SIZE} from "${projectDir}/modules/kallisto/kallisto_insert_size"
include {PIZZLY} from "${projectDir}/modules/pizzly/pizzly"

include {SQUID} from "${projectDir}/modules/squid/squid_call"
include {SQUID_ANNOTATE} from "${projectDir}/modules/squid/squid_annotate"

include {SAMTOOLS_VIEW as SAMTOOLS_VIEW_SQUID} from "${projectDir}/modules/samtools/samtools_view"

include {STAR_FUSION as STAR_FUSION} from "${projectDir}/modules/star-fusion/star-fusion"

include {FASTQC} from "${projectDir}/modules/fastqc/fastqc"
include {FUSION_REPORT} from "${projectDir}/modules/fusion_report/fusion_report"
include {MULTIQC} from "${projectDir}/modules/multiqc/multiqc"


// log params
param_log()

if (params.download_data && !params.csv_input) {
    exit 1, "Data download was specified with `--download_data`. However, no input CSV file was specified with `--csv_input`. This is an invalid parameter combination. `--download_data` requires a CSV manifest. See `--help` for information."
}

if (params.pdx && params.gen_org == 'mouse') {
    exit 1, "PDX analysis was specified with `--pdx`. `--gen_org` was set to: ${params.gen_org}. This is an invalid parameter combination. `--gen_org` must == 'human' for PDX analysis."
}

if (params.gen_org == 'mouse') {
    exit 1, "This pipeline currently only supports human data analysis."
}

if (params.read_type == 'SE') {
    exit 1, "This pipeline supports only paired end data."
}

// prepare reads channel
if (params.csv_input) {

    ch_input_sample = extract_csv(file(params.csv_input, checkIfExists: true))
    
    ch_input_sample.map{it -> [it[0], [it[2], it[3]]]}.set{read_ch}
    ch_input_sample.map{it -> [it[0], it[1]]}.set{meta_ch}
    
} else if (params.concat_lanes){
  
    read_ch = Channel
            .fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true, flat:true )
            .map { file, file1, file2 -> tuple(getLibraryId(file), file1, file2) }
            .groupTuple()
  
    // if channel is empty give error message and exit
    read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

} else {

    read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )

    // if channel is empty give error message and exit
    read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

}

// main workflow
workflow RNA_FUSION {

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

    // Step 0: Concat local Fastq files if required.
    if (params.concat_lanes && !params.csv_input){
        CONCATENATE_READS_PE(read_ch)
        read_ch = CONCATENATE_READS_PE.out.concat_fastq
    }

    GUNZIP(read_ch)

    FASTQ_PAIR(GUNZIP.out.gunzip_fastq)

    FASTQC(GUNZIP.out.gunzip_fastq)

    // Step 1a: Xenome if PDX data used.
    ch_XENOME_CLASSIFY_multiqc = Channel.empty() //optional log file. 
    if (params.pdx){
        // Xenome Classification
        XENOME_CLASSIFY(FASTQ_PAIR.out.paired_fastq)
        ch_XENOME_CLASSIFY_multiqc = XENOME_CLASSIFY.out.xenome_stats //set log file for multiqc

        // Xenome Read Sort
        FASTQ_SORT_HUMAN(XENOME_CLASSIFY.out.xenome_fastq, 'human')
        FASTQ_SORT_MOUSE(XENOME_CLASSIFY.out.xenome_mouse_fastq, 'mouse') 
        fusion_tool_input = FASTQ_SORT_HUMAN.out.sorted_fastq

    } else { 
        fusion_tool_input = FASTQ_PAIR.out.paired_fastq
    }

    // Step 3: Callers:
    // arriba
    STAR_ARRIBA(fusion_tool_input, params.arriba_star_args)
    SORT_ARRIBA(STAR_ARRIBA.out.bam, '')
    INDEX_ARRIBA(SORT_ARRIBA.out.sorted_bam)
    arriba_input = SORT_ARRIBA.out.sorted_bam.join(INDEX_ARRIBA.out.bai)
    ARRIBA(arriba_input)

    // fusioncatcher
    FUSIONCATCHER(fusion_tool_input)

    // jaffa
    JAFFA(fusion_tool_input)

    // pizzly
    KALLISTO_QUANT(fusion_tool_input)
    KALLISTO_INSERT_SIZE(KALLISTO_QUANT.out.kallisto_abundance)
    pizzly_input = KALLISTO_QUANT.out.kallisto_fusions.join(KALLISTO_INSERT_SIZE.out.kallisto_insert_size)
    PIZZLY(pizzly_input)

    // squid
    STAR_SQUID(fusion_tool_input, params.squid_star_args)
    SAMTOOLS_VIEW_SQUID(STAR_SQUID.out.sam, '-Sb', '_chimeric') // NOTE: The sam file from STAR_SQUID contains chimeric reads. Per STAR passed arguments. 
    SORT_SQUID(SAMTOOLS_VIEW_SQUID.out.bam, '')
    squid_input = STAR_SQUID.out.bam_sorted.join(SORT_SQUID.out.sorted_bam )
    SQUID(squid_input)
    SQUID_ANNOTATE(SQUID.out.squid_fusions)
    
    // star-fusion
    STAR_FUSION(fusion_tool_input)

    // Step 4: Fusion Reporter
    fusion_report_input = ARRIBA.out.arriba_fusions.join(FUSIONCATCHER.out.fusioncatcher_fusions).join(JAFFA.out.jaffa_fusions).join(PIZZLY.out.pizzly_fusions).join(SQUID_ANNOTATE.out.squid_fusions_annotated).join(STAR_FUSION.out.star_fusion_fusions)
    FUSION_REPORT(fusion_report_input)

    // Step 5: MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FUSION_REPORT.out.summary_fusions_mq.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.quality_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_XENOME_CLASSIFY_multiqc.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
}
