#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/emase.nf"
include {param_log} from "${projectDir}/bin/log/emase.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {BOWTIE} from "${projectDir}/modules/bowtie/bowtie"
include {SAMTOOLS_VIEW} from "${projectDir}/modules/samtools/samtools_view"
include {GBRS_BAM2EMASE} from "${projectDir}/modules/gbrs/gbrs_bam2emase"
include {GBRS_COMPRESS as GBRS_COMPRESS_SE;
         GBRS_COMPRESS as GBRS_COMPRESS_PE} from "${projectDir}/modules/gbrs/gbrs_compress"
include {GBRS_QUANTIFY} from "${projectDir}/modules/gbrs/gbrs_quantify"
include {EMASE_RUN} from "${projectDir}/modules/emase/emase_run"

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
        read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
    }
}

    // SE READS IN THIS CONTEXT NEEDS TESTING. 


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
            read_ch = CONCATENATE_READS_SE.out.concat_fastq
        }
    }

    // CONCAT READS IN THIS CONTEXT NEEDS TESTING. BOTH SE AND PE. 

    BOWTIE(read_ch)
    SAMTOOLS_VIEW(BOWTIE.out.sam, '-bS')
    GBRS_BAM2EMASE(SAMTOOLS_VIEW.out.bam)
    GBRS_COMPRESS_SE(GBRS_BAM2EMASE.out.emase_h5, '')

    if (params.read_type == 'PE'){
        gbrs_compress_pairedReads_input = GBRS_COMPRESS_SE.out.compressed_emase_h5
                                            .groupTuple()
                                            .map { sampleID, reads -> tuple( sampleID, reads.sort{it.name} ) }
        // collect GBRS compression. by sample ID then map and sort tuple to [sampleID, [R1, R2]]
        // NOTE THIS NEEDS TO BE CHECKED FOR MULTI-SAMPLE MIXING. 
        GBRS_COMPRESS_PE(gbrs_compress_pairedReads_input, 'merged')
        gbrs_quantify_input = GBRS_COMPRESS_PE.out.compressed_emase_h5
    } else {
        gbrs_quantify_input = GBRS_COMPRESS_SE.out.compressed_emase_h5
    }

    GBRS_QUANTIFY(gbrs_quantify_input)
    // NOTE: gbrs quantify is a wrapper around the `run-emase` code.

    EMASE_RUN(gbrs_quantify_input)
    // CHECK OUTPUT FROM THIS AGAINST GBRS. IS IT THE SAME? 
    
    // gbrs_compress_pairedReads_input.view()

}