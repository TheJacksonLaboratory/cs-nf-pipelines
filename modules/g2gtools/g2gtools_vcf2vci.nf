process G2GTOOLS_VCF2VCI {
    tag "$strain"

    cpus 8
    memory 5.GB
    time '02:30:00'


    container 'quay.io/jaxcompsci/g2gtools:0.2.9'

    publishDir "${params.pubdir}/g2gtools", pattern: '*.vci.gz*', mode:'copy'
    publishDir "${params.pubdir}/g2gtools", pattern: '*.errors.vcf', mode:'copy'

    input:
    val(strain)

    output:
    tuple val(strain), path("*.vci.gz"), path("*vci.gz.tbi"), emit: vci_tbi
    path("*.errors.vcf"), optional: true

    script:

    run_diploid = params.diploid ? '--diploid' : ''
    run_keep = params.keep_fails ? '--keep' : ''
    pass_only = params.pass_only ? '--pass' : ''
    debug_run = params.debug ? '--debug' : ''

    """
    /g2gtools/bin/g2gtools vcf2vci -p ${task.cpus} ${run_diploid} ${run_keep} ${pass_only} ${params.quality_filter} ${debug_run} -i ${params.snp_vcf} -i ${params.indel_vcf} -f ${params.primary_reference_fasta} -s ${strain} -o ${strain}.${params.genome_version}.vci
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.vci.gz
    touch ${strain}.${params.genome_version}.vci.gz.tbi
    """

}


/*
    Create VCI file from VCF file(s)

    Usage: g2gtools vcf2vci [-i <indel VCF file>]* -s <strain> -o <output file> [options]

    Required Parameters:
        -i, --vcf <vcf_file>             VCF file name
        -f, --fasta <Fasta File>         Fasta file matching VCF information
        -s, --strain <Strain>            Name of strain (column in VCF file)
        -o, --output <Output file>       VCI file name to create

    Optional Parameters:
        -p, --num-processes <number>     The number of processes to use, defaults to the number of cores
        --diploid                        Create diploid VCI
        --keep                           Keep track of VCF lines that could not be converted to VCI file
        --pass                           Use only VCF lines that have a PASS for the filter value
        --quality                        Filter on quality, FI=PASS
        --no-bgzip                       DO NOT compress and index output

    Help Parameters:
        -h, --help                       Print the help and exit
        -d, --debug                      Turn debugging mode on (list multiple times for more messages)
*/