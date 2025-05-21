#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {ARIA_DOWNLOAD} from "${projectDir}/modules/utility_modules/aria_download"
include {CONCATENATE_READS_SAMPLESHEET} from "${projectDir}/modules/utility_modules/concatenate_reads_sampleSheet"

workflow FILE_DOWNLOAD {

    take:
        ch_input_sample

    main:

    /* 
        General note: 

        Input tuple expected from the CSV sheet: 
            it[0] is sample ID. 
            it[1] is metadata information. meta includes: [sampleID:'testSample_42', lane:'lane1', replicate:'NA', id:'testSample_42', size:1]. This comes from `extract_csv.nf`
            it[2] and it[3] are R1 and R2 if PE. it[3] is empty if SE. 

        All steps expect that sampleID is in position [0] of tuples. 
        
        merge_replicates is used in the ATAC workflow and is used to merge replicates of one sample.

    */

        if (params.read_type == 'PE') {
            aria_download_input = ch_input_sample
            .multiMap { it ->
                if (params.merge_replicates) {
                    sampleID   = it[1].sampleID+'_'+it[1].replicate
                } else {
                    sampleID   = it[1].sampleID
                    ind   = it[1].ind
                    sex   = it[1].sex
                }
                R1: tuple(sampleID, it[1], 'R1', it[2])
                R2: tuple(sampleID, it[1], 'R2', it[3])
            }
            .mix()
            group_size = 2
        } else {
            aria_download_input = ch_input_sample
            .multiMap { it ->
                if (params.merge_replicates) {
                    sampleID   = it[1].sampleID+'_'+it[1].replicate
                } else {
                    sampleID   = it[1].sampleID
                    ind   = it[1].ind
                    sex   = it[1].sex
                }
                R1: tuple(sampleID, it[1], 'R1', it[2])
            }
            group_size = 1
        }
        /* 
            remap the data to individual R1 / R2 tuples. 
            These individual tuples are then mixed to pass individual files to the downloader. 
            R1 vs. R2 is maintained in the mix. Order is irrelavent here as data are grouped
            by sampleID downstream. 
        */

        aria_download_input.view()

        // Download files. 
        ARIA_DOWNLOAD(aria_download_input)

        concat_input = ARIA_DOWNLOAD.out.file
                            .map { it ->
                                def meta = [:]
                                meta.sampleID   = it[1].sampleID
                                meta.ind        = it[1].ind
                                meta.sex        = it[1].sex
                                [it[0], it[1].lane, meta, it[2], it[3], it[1].size] // sampleID, laneID, meta, read_ID:[R1|R2], file, number_of_lanes
                            }
                            .map {  sampleID, laneID, meta, readID, file, size -> tuple( groupKey([sampleID, meta, readID], size), laneID, file )  }
                            .groupTuple() // controlled by group key: [sampleID, meta, read_ID] 
                            .map{ it -> tuple(it[0][0], it[1].size(), it[0][1], it[0][2], it[2])} // sampleID, num_lanes, meta, read_ID:[R1|R2], file
                            .branch{
                                concat: it[1] > 1
                                pass:  it[1] == 1
                            }
        /* 
            remap the downloaded files to exclude lane from meta, and group on sampleID, meta, and read_ID: R1|R2.
            The number of lanes in the grouped data is used to determine if concatenation is needed. 
            The branch statement makes a 'concat' set for concatenation and a 'pass' set that isn't concatenated.
            The branch is using it[1].size() from the preceding step, i.e., the list size of lanes for the sample.  

            Metadata inclusion here is for future expansion. As implimented above, metadata is redundant to sampleID in `it[0]`. 
            However, if additional metadata are added to sample sheets, those metadata can be added and tracked above.

            groupTuple size is dynamically defined by metadata field 'size' i.e., the number of lanes per sample. 

            See: https://www.nextflow.io/docs/latest/operator.html#grouptuple and the note about dynamic group size. 
    
        */


        no_concat_samples = concat_input.pass
                            .map{it -> tuple(it[0], it[1], it[2], it[3], it[4][0])} // sampleID, num_lanes, meta, read_ID:[R1|R2], file
        /* 
            this delists the the file in `it[4]` as it is a single fastq sample (i.e., non-concat samples).
        */

        // Concatenate samples as needed. 
        CONCATENATE_READS_SAMPLESHEET(concat_input.concat)

        read_meta_ch = CONCATENATE_READS_SAMPLESHEET.out.concat_fastq
        .mix(no_concat_samples)
        .groupTuple(by: [0,2], size: group_size) // sampleID, meta
        .map{it -> tuple(it[0], it[2], it[4].toSorted( { a, b -> a.getName() <=> b.getName() } ) ) }

        read_meta_ch.view()

        /*
            Mix concatenation files, with non-concat files. 'mix' allows for, all, some, or no files to have 
            gone through concatenation. 

            Reads are remapped to read_ch and meta is placed in meta_ch. Input tuples for existing modules 
            do not expect 'meta' in the tuple. Example expected input tuple: [sampleID, [reads]]
        */
    emit:
        read_meta_ch
}