#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================
                     mmrSVD pipeline
========================================================================

The Mouse Mutant Resource Structural Variant Detection (mmrSVD) pipeline,
developed by Brian Sanderson (brian.sanderson@jax.org) and 
Mike Lloyd (mike.lloyd@jax.org) for The Jackson Laboratory.

*/

if (params.workflow == "pacbio") {
	include {PACBIO} from "./workflows/pacbio"
}

if (params.workflow == "illumina") {
    include {ILLUMINA} from "./workflows/illumina"
}

if (params.workflow == "ont") {
    include {ONT} from "./workflows/ont"
}

workflow {
    if (params.workflow == "pacbio") {
        PACBIO()
    }

    if (params.workflow == "illumina") {
        ILLUMINA()
    }

    if (params.workflow == "ont") {
        ONT()
    }    
}