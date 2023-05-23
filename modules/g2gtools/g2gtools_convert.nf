process G2GTOOLS_CONVERT {
    tag "$strain"

    cpus 1
    memory 5.GB
    time '02:30:00'

    container 'quay.io/jaxcompsci/g2gtools:0.2.9'

    publishDir "${params.pubdir}/g2gtools", pattern: "*.${format}", mode:'copy'

    input:
    tuple val(strain), path(vci), path(tbi)
    path(input_file)
    val(format)
    val(reverse)

    output:
    tuple val(strain), path("*.${format}"), emit: coverted_file
    tuple val(strain), path("*.${format}.unmapped"), emit: unmapped_file

    script:

    debug_run = params.debug ? '--debug' : ''
    run_reverse = reverse ? '--reverse' : ''

    if (params.append_chromosomes.find() == true | params.append_chromosomes == null) {
        append_chroms = false
    } else {
        append_chroms = true
    }
    // params.append_chromosomes.find() looks at the string params.append_chromosomes.
    // If there is no string, it returns 'true'. Also, we check if the param is 'null' as no append is needed in that case. 

    """
    /g2gtools/bin/g2gtools convert ${debug_run} -i ${input_file} -c ${vci} --format ${format} ${run_reverse} -o ${strain}.${params.genome_version}.${format}

    if [ ${append_chroms} ]
    then
        chroms_to_append=\$(echo ${params.append_chromosomes} | tr "," "\\n")
        
        for chrom in \$chroms_to_append
        do
            grep -P "\${chrom}\\t" ${strain}.${params.genome_version}.${format}.unmapped >> ${strain}.${params.genome_version}.${format}
        done
        
    fi

    ### CHECK IF VCI HAS ENTRIES FOR ALL CHROMS
    ### IF CHROM IS MISSING FROM VCI CONVERSION CALLS
    ### ADD BACK IN ANYTHING FROM UNMAPPED THAT IS NOT PRESENT IN THE VCI
    ### IF A SINGLE VARIANT IS THERE, HOW DOES IT WORK? DOES WHOLE CHROM GO? 
    
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.${format}
    """

}

/*

    ----- tool tip ----------
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
    ----------------------


    Note regarding 'append_chromosomes': Sanger in does not include 'Y' or 'MT' in the VCF file (Y for biological reasons, MT for seq depth reasons). 
                                         Having the option to map to Y and MT should be available. This extra block of code appends grep 
                                         matched strings: "<STRING>'\t'" bottom to the GTF file.



*/