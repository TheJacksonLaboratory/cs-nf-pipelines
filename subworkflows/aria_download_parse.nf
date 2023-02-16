#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {ARIA_DOWNLOAD} from "${projectDir}/modules/utility_modules/aria_download"
include {CONCATENATE_READS_SAMPLESHEET} from "${projectDir}/modules/utility_modules/concatenate_reads_sampleSheet"

workflow FILE_DOWNLOAD {

    take:
        ch_input_sample

    main:
        if (params.read_type == 'PE') {
            aria_download_input = ch_input_sample
            .multiMap { it ->
                R1: tuple(it[0], it[1], 'R1', it[2])
                R2: tuple(it[0], it[1], 'R2', it[3])
            }
            .mix()
        } else {
            aria_download_input = ch_input_sample
            .multiMap { it ->
                R1: tuple(it[0], it[1], 'R1', it[2])
            }
            .mix()
        }
        /* 
            remap the data to individual R1 / R2 tuples. 
            These individual tuples are then mixed to pass individual files to the downloader. 
            R1 vs. R2 is maintained in the mix. Order is irrelavent here as data are grouped
            by sampleID downstream. 
        */

        // Download files. 
        ARIA_DOWNLOAD(aria_download_input)

        concat_input = ARIA_DOWNLOAD.out.file
                            .map { it ->
                                def meta = [:]
                                meta.sample   = it[1].sample
                                meta.patient  = it[1].patient
                                meta.sex      = it[1].sample
                                meta.status   = it[1].sample
                                meta.id       = it[1].id

                                [it[0], it[1].lane, meta, it[2], it[3]]
                            }    
                            .groupTuple(by: [0,2,3])
                            .map{ it -> tuple(it[0], it[1].size(), it[2], it[3], it[4])}
                            .branch{
                                concat: it[1] > 1
                                pass:  it[1] == 1
                            }
        /* 
            remap the downloaded files to exclude lane from meta, and group on sampleID, meta, and read_num: R1|R2.
            The number of lanes in the grouped data is used to determine if concatenation is needed. 
            The branch statement makes a 'concat' set for concatenation and a 'pass' set that isn't concatenated. 
        */

        no_concat_samples = concat_input.pass
                            .map{it -> tuple(it[0], it[1], it[2], it[3], it[4][0])}
        /* 
            this delists the single fastq samples (i.e., non-concat samples).
        */

        // Concatenate samples as needed. 
        CONCATENATE_READS_SAMPLESHEET(concat_input.concat)

        read_meta_ch = CONCATENATE_READS_SAMPLESHEET.out.concat_fastq
        .mix(no_concat_samples)
        .groupTuple(by: [0,2])
        .map{it -> tuple(it[0], it[2], it[4].toSorted( { a, b -> a.getName() <=> b.getName() } ) ) }

        /*
            Mix concatenation files, with non-concat files. 'mix' allows for, all, some, or no files to have 
            gone through concatenation. 

            Reads are remapped to read_ch and meta is placed in meta_ch. Input tuples for existing modules 
            do not expect 'meta' in the tuple. Example expected input tuple: [sampleID, [reads]]
        */
    emit:
        read_meta_ch
}