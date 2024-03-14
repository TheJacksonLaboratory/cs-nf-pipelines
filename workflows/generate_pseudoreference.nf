#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/generate_pseudoreference.nf"
include {param_log} from "${projectDir}/bin/log/generate_pseudoreference.nf"
include {FILTER_GTF} from "${projectDir}/modules/utility_modules/filter_gtf_biotypes"
include {G2GTOOLS_VCF2VCI} from "${projectDir}/modules/g2gtools/g2gtools_vcf2vci"
include {G2GTOOLS_PATCH} from "${projectDir}/modules/g2gtools/g2gtools_patch"
include {G2GTOOLS_TRANSFORM} from "${projectDir}/modules/g2gtools/g2gtools_transform"
include {SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_G2GTOOLS;
         SAMTOOLS_FAIDX as SAMTOOLS_FAIDX} from "${projectDir}/modules/samtools/samtools_faidx"
include {G2GTOOLS_CONVERT} from "${projectDir}/modules/g2gtools/g2gtools_convert"
include {APPEND_DROPPED_CHROMS} from "${projectDir}/modules/python/append_dropped_chroms"
include {G2GTOOLS_GTF2DB} from "${projectDir}/modules/g2gtools/g2gtools_gtf2db"
include {G2GTOOLS_EXTRACT as G2GTOOLS_EXTRACT_GENES;
         G2GTOOLS_EXTRACT as G2GTOOLS_EXTRACT_TRANSCRIPTS;
         G2GTOOLS_EXTRACT as G2GTOOLS_EXTRACT_EXONS} from "${projectDir}/modules/g2gtools/g2gtools_extract"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

if (params.region != '' && params.bed != '') {
    exit 1, "Both REGION: ${params.region} and BED: ${params.bed} can not be set. Use only one of these options at a time."
}

// main workflow
workflow GENERATE_PSEUDOREFERENCE  {

    if (params.gtf_biotype_include) {
        FILTER_GTF()
        gtf_file = FILTER_GTF.out.filtered_gtf
    } else {
        gtf_file = params.primary_reference_gtf
    }

    Channel
    .of( params.strain.split(',') )
    .set { strain }
    /* 
       The above makes a channel by splitting a comma delim string and passes that into a forked channel. 
       In this case the string is a list of strains. 

       Note that if a simple list is passed as input, 
       Nextflow will keep a list together as a basic data type, 
       and not split into multiple processes. 

        FROM NEXTFLOW:  
        Workflow inputs are always channels by definition. 
        If a basic data type is provided instead, such as a number, string, list, etc, 
        it is implicitly converted to a value channel.

        See also: https://www.nextflow.io/docs/latest/channel.html#of
    */

    G2GTOOLS_VCF2VCI(strain)
    G2GTOOLS_PATCH(G2GTOOLS_VCF2VCI.out.vci_tbi)
    transform_input = G2GTOOLS_PATCH.out.patched_fasta.join(G2GTOOLS_VCF2VCI.out.vci_tbi)
    G2GTOOLS_TRANSFORM(transform_input)
    SAMTOOLS_FAIDX_G2GTOOLS(G2GTOOLS_TRANSFORM.out.final_fasta)
    G2GTOOLS_CONVERT(G2GTOOLS_VCF2VCI.out.vci_tbi, gtf_file, 'gtf', false)
    
    if (params.append_chromosomes) {
        APPEND_DROPPED_CHROMS(G2GTOOLS_VCF2VCI.out.vci_tbi.join(G2GTOOLS_CONVERT.out.unmapped_file).join(G2GTOOLS_CONVERT.out.coverted_file))
        gtf2db_input = APPEND_DROPPED_CHROMS.out.appended_gtf
    } else {
        gtf2db_input = G2GTOOLS_CONVERT.out.coverted_file
    }
    
    G2GTOOLS_GTF2DB(gtf2db_input)
    extract_input = G2GTOOLS_TRANSFORM.out.final_fasta.join(G2GTOOLS_GTF2DB.out.db)
    G2GTOOLS_EXTRACT_GENES(extract_input, 'genes')
    G2GTOOLS_EXTRACT_TRANSCRIPTS(extract_input, 'transcripts')
    G2GTOOLS_EXTRACT_EXONS(extract_input, 'exons')

    /*
    For each STRAIN the following steps were run: 
        1. Convert VCF to VCI (chain file equivalent)
        2. Path SNPs into reference. 
        3. Transform InDELs into patched reference. 
        5. Convert the reference GTF to strain specific GTF.
        6. Convert strain specific GTF to database format.
        7. Extract sequence from strain specific fasta: 
            a. genes.
            b. transcripts.
            c. exons.
    */
    
    faidx_input = ['primary_strain', params.primary_reference_fasta]

    SAMTOOLS_FAIDX(faidx_input)
    // GBRS requies 'ref.fa.idx' which is fasta index of the primary refernce. 
    // The index of that file is easily done here. 
    // SAMTOOLS_FAIDX requires an input tuple. It is dummied here to strain = 'primary_strain', fasta = primary_reference_fasta
}