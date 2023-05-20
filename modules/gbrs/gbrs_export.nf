process GBRS_EXPORT {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'

    container 'quay.io/jaxcompsci/gbrs_py3:feature_py3-b362dec'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/gbrs' : 'gbrs' }", pattern: "*.tsv", mode: 'copy'

    input:
    tuple val(sampleID), path(interpolated_genoprobs)

    output:
    tuple val(sampleID), file("*.gbrs.interpolated.genoprobs.tsv"), emit: interpolated_genoprobs_tsv


    script:

    """
    gbrs export \
        -i ${interpolated_genoprobs} \
        -s ${params.gbrs_strain_list} \
        -g ${params.genotype_grid} \
        -o ${sampleID}.gbrs.interpolated.genoprobs.tsv
    """

    stub:
    """
    touch ${sampleID}.gbrs.interpolated.genoprobs.tsv
    """
}

/*
 Usage: gbrs export [OPTIONS]

 export to GBRS quant format

╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --genoprob-file  -i      FILE     genotype probability file [default: None] [required]                                                                                                                                                                          │
│ *  --strains        -s      TEXT     strain, either one per -s option, i.e. -h A -h B -h C, or a shortcut -s A,B,C [default: None] [required]                                                                                                                      │
│    --grid-file      -g      FILE     grid file (i.e, ref.genome_grid.64k.txt) [default: None]                                                                                                                                                                      │
│    --output         -o      FILE     output file in GBRS quant format [default: None]                                                                                                                                                                              │
│    --verbose        -v      INTEGER  specify multiple times for more verbose output [default: 0]                                                                                                                                                                   │
│    --help                            Show this message and exit.                                                                                                                                                                                                   │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯


    // From the script, default param setttings: 
    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-i",
            "--genoprob-file",
            action="store",
            dest="gprobfile",
            type=str,
            required=True,
        )
        parser.add_argument(
            "-s", "--strains", action="store", dest="strains", type=str, required=True
        )
        parser.add_argument(
            "-g", "--grid-file", action="store", dest="gridfile", type=str, default=None
        )
        parser.add_argument(
            "-f", "--in-numpy-format", action="store_true"  # file format
        )
        parser.add_argument(
            "-c", "--by-chromosomes", action="store_true"  # row dimension
        )
        parser.add_argument(
            "-p", "--by-strains", action="store_false"  # column dimension
        )
        parser.add_argument(
            "-o", action="store", dest="outfile", type=str, default=None
        )
        parser.add_argument(
            "-v", "--verbose", help="Toggle DEBUG verbosity", action="store_true"
        )
        return parser.parse_args()
*/
