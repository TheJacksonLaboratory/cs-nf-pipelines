process G2GTOOLS_CONVERT {

    cpus 1
    memory {15.GB * task.attempt}
    time {10.hour * task.attempt}
    errorStrategy 'retry' 
    maxRetries 1

    container 'quay.io/jaxcompsci/g2gtools:0.2.9'

    publishDir "${params.pubdir}/g2gtools", pattern: "*.${format}", mode:'copy'

    input:
    tuple val(strain), path(vci), path(tbi)
    path(input_file)
    val(format)
    val(reverse)

    output:
    tuple val(strain), path("*.${format}"), emit: coverted_file

    script:

    debug_run = params.debug ? '--debug' : ''
    run_reverse = reverse ? '--reverse' : ''

    """
    /g2gtools/bin/g2gtools convert ${debug_run} -i ${input_file} -c ${vci} --format ${format} ${run_reverse} -o ${strain}.${params.genome_version}.${format}
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.${format}
    """

}

/*

    Liftover coordinates of bam|sam|gtf|bed files

    Usage: g2gtools convert -c <VCI file> -i <input file> [options]

    Required Parameters:
        -i, --input <input file>         Input file to liftover
        -c, --vci <VCI file>             VCI file

    Optional Parameters:
        -f, --format <bam|sam|gtf|bed>   Input file format
        -o, --output <output file>       Name of output file
        --reverse                        Reverse the direction of the conversion (Custom -> Ref)

    Help Parameters:
        -h, --help                       Print the help and exit
        -d, --debug                      Turn debugging mode on (list multiple times for more messages)


    Note:
        If no output file is specified [-o], the converted information will be redirected
        to standard out.

*/