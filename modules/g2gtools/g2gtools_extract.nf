process G2GTOOLS_EXTRACT {

    cpus 1
    memory 6.GB
    time '02:30:00'


    container 'quay.io/jaxcompsci/g2gtools:0.2.9'

    publishDir "${params.pubdir}/g2gtools", pattern: '*.fa', mode:'copy'

    input:
    tuple val(strain), path(final_fasta), path(db)
    val(extract_type)

    output:
    tuple val(strain), path("*.fa"), emit: extracted_fasta

    script:

    debug_run = params.debug ? '--debug' : ''

    """
    /g2gtools/bin/g2gtools extract ${debug_run} -i ${final_fasta} -db ${db} --${extract_type} > ${strain}.${params.genome_version}.${extract_type}.fa
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.${extract_type}.fa
    """

}
/* 
NOTE: The above script is hard-coded for extraction of regions from database files. Additional options are available. 
      At this time, only db to fasta is used. If additional functionality is required from this module, additional coding, 
      and variables will be required. See help text below for these options / input formats. Additional conditional logic 
      or sub-modules of the extract function may also be required to avoid parameter conflicts. 
*/ 

/*
    Extract subsequence from a fasta file given a region

    Usage: g2gtools extract -i <Fasta file> [-r <region> | -b <BED file> | -db <Database> | -id <seqid>] [options]

    Required Parameters:
        -i, --input <Fasta file>         Fasta file to extract subsequence from, BGZIP files supported
        -b, --bed <BED file>             BED file -OR-
        -id, --identifier <identifier>   Fasta identifier -OR-
        -r, --region <seqid:start-end>   Region to extract -OR-
        -db, --database <DB file>        Database file
             --transcripts                 For use with -db, extracts transcripts (default)
             --exons                       For use with -db, extracts exons
             --genes                       For use with -db, extracts genes

    Optional Parameters:
        -c, --vci <VCI file>             Input VCI file, matching input Fasta file
        -R, --vci-reverse                Reverse the direction of liftover
        --complement                     Complement the extracted sequence
        --reverse                        Reverse the extracted sequence
        --reverse-complement             Reverse complement the extracted sequence
        --raw                            Just shows the extracted sequence

    Help Parameters:
        -h, --help                       Print the help and exit
        -d, --debug                      Turn debugging mode on (list multiple times for more messages)

    Note:
        Locations specified on command line are 1-based coordinates.
        Locations specified via BED file are 0-based coordinates.
*/