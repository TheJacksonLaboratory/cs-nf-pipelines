process G2GTOOLS_VCF2VCI {
    tag "$strain"

    cpus 8
    memory 5.GB
    time '02:30:00'


    container 'quay.io/jaxcompsci/g2gtools:74926ad'

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
    g2gtools vcf2vci -p ${task.cpus} ${run_diploid} ${run_keep} ${pass_only} ${params.quality_filter} ${debug_run} -i ${params.snp_vcf} -i ${params.indel_vcf} -f ${params.primary_reference_fasta} -s ${strain} -o ${strain}.${params.genome_version}.vci
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.vci.gz
    touch ${strain}.${params.genome_version}.vci.gz.tbi
    """

}

/*
 Usage: g2gtools vcf2vci [OPTIONS]

 Create VCI file from VCF file(s)

╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --vcf       -i      FILE     VCF files can seperate files by "," or have multiple -i [default: None] [required]                                                                                                                                                                    │
│ *  --fasta     -f      FILE     Fasta file matching VCF information [default: None] [required]                                                                                                                                                                                        │
│ *  --strain    -s      TEXT     Name of strain/sample (column in VCF file) [default: None] [required]                                                                                                                                                                                 │
│    --output    -o      FILE     Name of output file [default: None]                                                                                                                                                                                                                   │
│    --diploid   -d               Create diploid VCI file                                                                                                                                                                                                                               │
│    --keep      -k               Keep track of VCF lines that could not be converted to VCI file                                                                                                                                                                                       │
│    --pass      -p               Use only VCF lines that have a PASS for the filter value                                                                                                                                                                                              │
│    --quality   -q               Filter on quality, FI=PASS                                                                                                                                                                                                                            │
│    --no-bgzip  -z               DO NOT compress and index output                                                                                                                                                                                                                      │
│    --verbose   -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                                           │
│    --help                       Show this message and exit.                                                                                                                                                                                                                           │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
*/
