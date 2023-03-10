process G2GTOOLS_PATCH {

    cpus 8
    memory 8.GB
    time '02:30:00'

    container 'quay.io/jaxcompsci/g2gtools:0.2.9'

    publishDir "${params.pubdir}/g2gtools", pattern: '*.patched.fa', mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(strain), path(vci), path(tbi)

    output:
    tuple val(strain), path("*.patched.fa"), emit: patched_fasta

    script:

    debug_run = params.debug ? '--debug' : ''

    """
    /g2gtools/bin/g2gtools patch -p ${task.cpus} ${params.region} ${params.bed} ${debug_run} -i ${params.primary_reference_fasta} -c ${vci} -o ${strain}.${params.genome_version}.patched.fa
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.patched.fa
    """

}

/*
    Patch SNPs onto the reference sequence

    Usage: g2gtools patch -i <Fasta file> -c <VCI file> [options]

    Required Parameters:
        -i, --input <Fasta file>         Reference fasta file to extract and patch on
        -c, --vci <VCI file>             VCI file

    Optional Parameters:
        -p, --num-processes <number>     The number of processes to use, defaults to the number of cores
        -o, --output <Output file>       Name of output fasta file that contains SNP-patched sequences
        -r, --region <seqid:start-end>   Region to extract (cannot be used with -b)
        -b, --bed <BED file>             BED file (cannot be used with -r)
        --bgzip                          Compress and index output (can only be used with -o)

    Help Parameters:
        -h, --help                       Print the help and exit
        -d, --debug                      Turn debugging mode on (list multiple times for more messages)

    Note:
        --bgzip can be potentially slow depending on the size of the file
*/