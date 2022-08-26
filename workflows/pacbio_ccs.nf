#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules

include {BUILDPBMM2INDEX} from "../modules/pbmm2"

include {PBMM2MAPCCS} from "../modules/pbmm2"

include {DISCOVERTANDEM} from "../modules/pbsv"

include {CALLCCS} from "../modules/pbsv"

workflow PACBIO_CCS {
    // Step 1: Prepare index
    BUILDPBMM2INDEX(params.names, params.fasta)

    // Step 2: Map CCS reads to indexed genome
    PBMM2MAPCCS(params.names, params.fastq1, BUILDPBMM2INDEX.out.pbmm2_index)

    // Step 3: Call SV with PBSV
    DISCOVERTANDEM(PBMM2MAPCCS.out.pbmm2_ccs)
    CALLCCS(DISCOVERTANDEM.out.pbsv_svsig, params.fasta)

    // Step 4: Call SV with Sniffles
    SNIFFLES(PBMM2MAPCCS.out.pbmm2_ccs)
}