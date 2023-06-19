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

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow -c /path/to/params.config run /path/to/mmrSVD/main.nf -profile slurm,singularity --genome mm10

    """
}

// Show help message
if (params.help) exit 0, helpMessage()

params.fasta = params.genome ? params.genomes[params.genome].fasta ?: null : null

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