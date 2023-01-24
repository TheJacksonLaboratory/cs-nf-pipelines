process G2GTOOLS_TRANSFORM {

    cpus 8
    memory 25.GB
    time '01:30:00'

    container 'quay.io/jaxcompsci/g2gtools:0.2.9'

    publishDir "${params.pubdir}/g2gtools", pattern: '*.fa', mode:'copy'

    input:
    tuple val(strain), path(patched_fasta), path(vci), path(tbi)

    output:
    tuple val(strain), path("*.fa"), emit: final_fasta

    script:

    debug_run = params.debug ? '--debug' : ''

    """
    /g2gtools/bin/g2gtools transform -p ${task.cpus} ${params.region} ${params.bed} ${debug_run} -i ${patched_fasta} -c ${vci} -o ${strain}.${params.genome_version}.fa
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.fa
    """

}

/*
    Incorporate indels onto the input sequence

    Usage: g2gtools transform -i <Fasta file> -c <VCI file> [options]

    Required Parameters:
        -i, --input <Fasta file>         Reference fasta file to extract and transform
        -c, --vci <VCI file>             VCI file

    Optional Parameters:
        -p, --num-processes <number>     The number of processes to use, defaults to the number of cores
        -o, --output <Output file>       Name of output fasta file
        -r, --region <seqid:start-end>   Region to extract (cannot be used with -b)
        -b, --bed <BED file>             BED file (cannot be used with -r)
        --bgzip                          Compress and index output fasta file (can only be used with -o)

    Help Parameters:
        -h, --help                       Print the help and exit
        -d, --debug                      Turn debugging mode on (list multiple times for more messages)

    Note:
        --bgzip can be potentially slow depending on the size of the file
*/