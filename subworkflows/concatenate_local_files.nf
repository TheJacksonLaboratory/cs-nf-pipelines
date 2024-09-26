#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {CONCATENATE_READS_SAMPLESHEET} from "${projectDir}/modules/utility_modules/concatenate_reads_sampleSheet"

workflow CONCATENATE_LOCAL_FILES {

    take:
        ch_input_sample

    main:

        if (params.read_type == 'PE') {
            temp_map = ch_input_sample
            .multiMap { it ->
                def meta = [:]
                if (params.merge_replicates) {
                    meta.sampleID   = it[1].replicate != 'NA' ? it[1].sampleID+'_'+it[1].replicate : it[1].sampleID
                    meta.replicate  = it[1].replicate
                    meta.baseSampleID = it[1].sampleID
                } else {
                    meta.sampleID   = it[1].sampleID
                }
                R1: tuple(meta.sampleID, meta.lane, meta, 'R1', it[2])
                R2: tuple(meta.sampleID, meta.lane, meta, 'R2', it[3])
            }
            .mix()
            .groupTuple(by: [0,2,3])
            .map{ it -> tuple(it[0], it[1].size(), it[2], it[3], it[4]) } // sampleID, num_lanes, meta, read_ID:[R1|R2], file
            
            concat_input = temp_map
            .branch {
                concat: it[1] > 1
                pass:  it[1] == 1
            }
            group_size = 2
        } else {

            temp_map = ch_input_sample
            .multiMap { it ->
                def meta = [:]
                if (params.merge_replicates) {
                    meta.sampleID   = it[1].replicate != 'NA' ? it[1].sampleID+'_'+it[1].replicate : it[1].sampleID
                    meta.replicate  = it[1].replicate
                    meta.baseSampleID = it[1].sampleID
                } else {
                    meta.sampleID   = it[1].sampleID
                }
                R1: tuple(meta.sampleID, meta.lane, meta, 'R1', it[2])
            }
            .groupTuple(by: [0,2,3])
            .map{ it -> tuple(it[0], it[1].size(), it[2], it[3], it[4]) } // sampleID, num_lanes, meta, read_ID:[R1], file
            
            concat_input = temp_map
            .branch {
                concat: it[1] > 1
                pass:  it[1] == 1
            }
            group_size = 1 
        }

        no_concat_samples = concat_input.pass
                            .map{it -> tuple(it[0], it[1], it[2], it[3], it[4][0])} // sampleID, num_lanes, meta, read_ID:[R1|R2], file

        /* 
            this delists the the file in `it[4]` as it is a single fastq sample (i.e., non-concat samples).

        */

        CONCATENATE_READS_SAMPLESHEET(concat_input.concat)

        read_meta_ch = CONCATENATE_READS_SAMPLESHEET.out.concat_fastq
        .mix(no_concat_samples)
        .groupTuple(by: [0,2], size: group_size) // sampleID, meta
        .map{it -> tuple(it[0], it[2], it[4].toSorted( { a, b -> file(a).getName() <=> file(b).getName() } ) ) }
       
        /*
            Mix concatenation files, with non-concat files. 'mix' allows for, all, some, or no files to have 
            gone through concatenation. 

            Reads are remapped to read_ch and meta is placed in meta_ch. Input tuples for existing modules 
            do not expect 'meta' in the tuple. Example expected input tuple: [sampleID, [reads]]
        */

    emit:
        read_meta_ch

}