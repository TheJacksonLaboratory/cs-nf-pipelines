#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/emase.nf"
include {param_log} from "${projectDir}/bin/log/emase.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId_emase.nf"
include {extract_csv} from "${projectDir}/bin/shared/extract_csv.nf"
include {FILE_DOWNLOAD} from "${projectDir}/subworkflows/aria_download_parse"
include {CONCATENATE_LOCAL_FILES} from "${projectDir}/subworkflows/concatenate_local_files"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {RUN_EMASE} from "${projectDir}/subworkflows/run_emase"

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
        ch_input_sample.map{it -> [it[0], it[2], 'R1']}.set{r1}
        ch_input_sample.map{it -> [it[0], it[3], 'R2']}.set{r2}
        read_ch = r1.mix(r2)
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
        temp_read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
        temp_read_ch.map{it -> [it[0], it[1][0], 'R1']}.set{r1}
        temp_read_ch.map{it -> [it[0], it[1][1], 'R2']}.set{r2}
        read_ch = r1.mix(r2)
    }
    else if (params.read_type == 'SE'){
        temp_read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
        temp_read_ch.map{it -> [it[0], it[1][0], 'R1']}.set{read_ch}
    }
    // if channel is empty give error message and exit
    read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern} and file extension: ${params.extension}"}
}

// main workflow
workflow EMASE {
    // Step 0: Download data and concat Fastq files if needed. 
    if (params.download_data){
        FILE_DOWNLOAD(ch_input_sample)

        if (params.read_type == 'PE'){
            FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[2][0], 'R1']}.set{r1}
            FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[2][1], 'R2']}.set{r2}
            read_ch = r1.mix(r2)
        } else if (params.read_type == 'SE'){
            FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[2][0], 'R1']}.set{read_ch}
        }

        FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

    // Step 00: Concat local Fastq files from CSV input if required.
    if (!params.download_data && params.csv_input){
        CONCATENATE_LOCAL_FILES(ch_input_sample)
        
        if (params.read_type == 'PE'){
            CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[2][0], 'R1']}.set{r1}
            CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[2][1], 'R2']}.set{r2}
            read_ch = r1.mix(r2)
        } else if (params.read_type == 'SE'){
            CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[2][0], 'R1']}.set{read_ch}
        }

        CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

    // Step 00: Concatenate Fastq files if required. 
    if (params.concat_lanes && !params.csv_input){
        if (params.read_type == 'PE'){
            CONCATENATE_READS_PE(read_ch)
            temp_read_ch = CONCATENATE_READS_PE.out.concat_fastq
            temp_read_ch.map{it -> [it[0], it[1][0], 'R1']}.set{r1}
            temp_read_ch.map{it -> [it[0], it[1][1], 'R2']}.set{r2}
            read_ch = r1.mix(r2)
        } else if (params.read_type == 'SE'){
            CONCATENATE_READS_SE(read_ch)
            temp_read_ch = CONCATENATE_READS_SE.out.concat_fastq
            temp_read_ch.map{it -> [it[0], it[1], 'R1']}.set{read_ch}
        }
    }

    RUN_EMASE(read_ch)
    // workflow found in: subworkflows/run-emase
    // workflow run as subworkflow due to re-use in GBRS workflow. 

}
