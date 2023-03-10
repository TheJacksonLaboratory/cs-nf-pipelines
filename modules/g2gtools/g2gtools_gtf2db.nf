process G2GTOOLS_GTF2DB {

    cpus 1
    memory 1.GB
    time '02:30:00'


    container 'quay.io/jaxcompsci/g2gtools:0.2.9'

    publishDir "${params.pubdir}/g2gtools", pattern: '*.gtf.db', mode:'copy'

    input:
    tuple val(strain), path(gtf)

    output:
    tuple val(strain), path("*.gtf.db"), emit: db

    script:

    debug_run = params.debug ? '--debug' : ''

    """
    /g2gtools/bin/g2gtools gtf2db ${debug_run} -i ${gtf} -o ${strain}.${params.genome_version}.gtf.db
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.gtf.db
    """

}

/*
    Convert a GTF file to a G2G DB file

    Usage: g2gtools gtf2db -i <GTF file> -o <G2G DB file> [options]

    Required Parameters:
        -i, --input <GTF file>           GTF file to parse
        -o, --output <G2G DB file>       G2G DB file to create

    Help Parameters:
        -h, --help                       Print the help and exit
        -d, --debug                      Turn debugging mode on (list multiple times for more messages)


*/