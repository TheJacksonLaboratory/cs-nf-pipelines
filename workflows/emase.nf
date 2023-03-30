#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/emase.nf"
include {param_log} from "${projectDir}/bin/log/emase.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {RUN_EMASE} from "${projectDir}/subworkflows/run-emase"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

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
        temp_read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
        temp_read_ch.map{it -> [it[0], it[1][0], 'R1']}.set{r1}
        temp_read_ch.map{it -> [it[0], it[1][1], 'R2']}.set{r2}
        read_ch = r1.mix(r2)
    }
    else if (params.read_type == 'SE'){
        temp_read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
        temp_read_ch.map{it -> [it[0], it[1][0], 'R1']}.set{read_ch}
    }
}


// if channel is empty give error message and exit
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

// main workflow
workflow EMASE {
    // Step 0: Concatenate Fastq files if required. 
    if (params.concat_lanes){
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
