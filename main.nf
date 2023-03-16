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

// fasta can be either given as a genome name in iGenomes or as a fasta file
// params.fasta will be changed only if it's not previously defined
/*
params.fasta = params.genome ? params.genomes[params.genome].fasta ?: null : null
ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : "null"
ch_fastq1 = params.fastq1 ? Channel.value(file(params.fastq1)) : null
ch_fastq2 = params.fastq2 ? Channel.value(file(params.fastq2)) : null
params.surv_dist = 1000
params.surv_supp = 1
params.surv_type = 1
params.surv_strand = 1
params.surv_min = 30
def sample_name = params.names
def abs_outdir = params.outdir
params.threads = 8
params.keep_intermediate = false
*/

params.fasta = params.genome ? params.genomes[params.genome].fasta ?: null : null

if (params.workflow == "pacbio_ccs") {
	include {PACBIO_CCS} from "./workflows/pacbio_ccs"
}

if (params.workflow == "illumina") {
    include {ILLUMINA} from "./workflows/illumina"
}

workflow {
	if (params.workflow == "pacbio_ccs") {
		PACBIO_CCS()
	}

    if (params.workflow == "illumina") {
        ILLUMINA()
    }
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Closing Info ~ ~ ~ ~ ~ ~
/*
workflow.onComplete {
    wfEnd = [:]
    wfEnd['Completed at'] = workflow.complete
    wfEnd['Duration']     = workflow.duration
    wfEnd['Exit status']  = workflow.exitStatus
    wfEnd['Success']      = workflow.success
    if(!workflow.success){
        wfEnd['!!Execution failed'] = ''
        wfEnd['.    Error']   = workflow.errorMessage
        wfEnd['.    Report']  = workflow.errorReport
    }
    Summary.show(wfEnd)
}
*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
