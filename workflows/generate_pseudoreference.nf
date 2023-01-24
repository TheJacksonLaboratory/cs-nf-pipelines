#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/generate_pseudoreference.nf"
include {param_log} from "${projectDir}/bin/log/generate_pseudoreference.nf"
include {G2GTOOLS_VCF2VCI} from "${projectDir}/modules/g2gtools/g2gtools_vcf2vci"
include {G2GTOOLS_PATCH} from "${projectDir}/modules/g2gtools/g2gtools_patch"
include {G2GTOOLS_TRANSFORM} from "${projectDir}/modules/g2gtools/g2gtools_transform"
include {G2GTOOLS_CONVERT} from "${projectDir}/modules/g2gtools/g2gtools_convert"
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
    // strain = .toSortedList()
    // strain

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
    G2GTOOLS_CONVERT(G2GTOOLS_VCF2VCI.out.vci_tbi, params.primary_reference_gtf, 'gtf', false)
    G2GTOOLS_GTF2DB(G2GTOOLS_CONVERT.out.coverted_file)
    extract_input = G2GTOOLS_TRANSFORM.out.final_fasta.join(G2GTOOLS_GTF2DB.out.db)
    G2GTOOLS_EXTRACT_GENES(extract_input, 'genes')
    G2GTOOLS_EXTRACT_TRANSCRIPTS(extract_input, 'transcripts')
    G2GTOOLS_EXTRACT_EXONS(extract_input, 'exons')
}